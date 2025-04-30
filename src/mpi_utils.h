#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include <mpi.h>
#include "simulation.h"

// Structure to hold MPI-related data
typedef struct {
    int rank;           // Process rank
    int size;           // Total number of processes
    int local_nx;       // Local grid size in x
    int local_ny;       // Local grid size in y
    int start_x;        // Starting x index for this process
    int start_y;        // Starting y index for this process
    int *send_counts;   // Number of elements to send to each process
    int *displs;        // Displacement for each process in gather/scatter
} MPI_Data;

// Function declarations
MPI_Data* initialize_mpi_data(int NX, int NY);
void free_mpi_data(MPI_Data *mpi_data);
void distribute_grid(Simulation *sim, MPI_Data *mpi_data);
void gather_grid(Simulation *sim, MPI_Data *mpi_data);
void exchange_ghost_cells(Simulation *sim, MPI_Data *mpi_data);

#endif // MPI_UTILS_H 