#include <stdlib.h>
#include <string.h>
#include "mpi_utils.h"

MPI_Data* initialize_mpi_data(int NX, int NY) {
    MPI_Data *mpi_data = malloc(sizeof(MPI_Data));
    if (!mpi_data) return NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_data->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_data->size);

    // Calculate local grid dimensions
    mpi_data->local_nx = NX;
    mpi_data->local_ny = NY / mpi_data->size;
    if (mpi_data->rank < NY % mpi_data->size) {
        mpi_data->local_ny++;
    }

    // Calculate starting indices
    mpi_data->start_y = 0;
    for (int i = 0; i < mpi_data->rank; i++) {
        int ny = NY / mpi_data->size;
        if (i < NY % mpi_data->size) ny++;
        mpi_data->start_y += ny;
    }
    mpi_data->start_x = 0;

    // Allocate arrays for gather/scatter operations
    mpi_data->send_counts = malloc(mpi_data->size * sizeof(int));
    mpi_data->displs = malloc(mpi_data->size * sizeof(int));

    // Calculate send counts and displacements
    int current_displ = 0;
    for (int i = 0; i < mpi_data->size; i++) {
        int ny = NY / mpi_data->size;
        if (i < NY % mpi_data->size) ny++;
        mpi_data->send_counts[i] = NX * ny;
        mpi_data->displs[i] = current_displ;
        current_displ += mpi_data->send_counts[i];
    }

    return mpi_data;
}

void free_mpi_data(MPI_Data *mpi_data) {
    if (!mpi_data) return;
    free(mpi_data->send_counts);
    free(mpi_data->displs);
    free(mpi_data);
}

void exchange_ghost_cells(Simulation *sim, MPI_Data *mpi_data) {
    int NX = sim->NX;
    int local_ny = mpi_data->local_ny;
    
    // Exchange ghost cells with neighbors
    if (mpi_data->rank > 0) {
        // Send bottom row to previous process
        MPI_Send(&sim->density[IX(0, 0, NX)], NX, MPI_DOUBLE, 
                 mpi_data->rank - 1, 0, MPI_COMM_WORLD);
        // Receive ghost cells from previous process
        MPI_Recv(&sim->density[IX(0, -1, NX)], NX, MPI_DOUBLE,
                 mpi_data->rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    if (mpi_data->rank < mpi_data->size - 1) {
        // Send top row to next process
        MPI_Send(&sim->density[IX(0, local_ny - 1, NX)], NX, MPI_DOUBLE,
                 mpi_data->rank + 1, 1, MPI_COMM_WORLD);
        // Receive ghost cells from next process
        MPI_Recv(&sim->density[IX(0, local_ny, NX)], NX, MPI_DOUBLE,
                 mpi_data->rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void distribute_grid(Simulation *sim, MPI_Data *mpi_data) {
    // Scatter the grid data to all processes
    MPI_Scatterv(sim->density, mpi_data->send_counts, mpi_data->displs,
                 MPI_DOUBLE, &sim->density[IX(0, 0, sim->NX)],
                 mpi_data->local_nx * mpi_data->local_ny,
                 MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void gather_grid(Simulation *sim, MPI_Data *mpi_data) {
    // Gather the grid data from all processes
    MPI_Gatherv(&sim->density[IX(0, 0, sim->NX)],
                mpi_data->local_nx * mpi_data->local_ny,
                MPI_DOUBLE, sim->density,
                mpi_data->send_counts, mpi_data->displs,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
} 