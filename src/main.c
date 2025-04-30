#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "simulation.h"
#include "grid.h"
#include "visualization.h"
#include "mpi_utils.h"

int main(int argc, char **argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Simulation grid parameters.
    int NX = 100;           // Number of grid cells in x-direction
    int NY = 100;           // Number of grid cells in y-direction
    double Lx = 1.0;        // Domain length in X
    double Ly = 1.0;        // Domain length in Y

    // Time-related parameters.
    double dt = 0.005;       // Time step size
    double T = 5.0;        // Total simulation time

    // Output frequency: output every 'output_interval' simulation steps.
    int output_interval = 1;

    // Initialize MPI data structure
    MPI_Data *mpi_data = initialize_mpi_data(NX, NY);
    if (!mpi_data) {
        fprintf(stderr, "Error initializing MPI data\n");
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // Create the simulation grid (only on rank 0)
    Grid *grid = NULL;
    if (mpi_data->rank == 0) {
        grid = create_grid(NX, NY, Lx, Ly);
        if (!grid) {
            fprintf(stderr, "Error creating grid\n");
            free_mpi_data(mpi_data);
            MPI_Finalize();
            return EXIT_FAILURE;
        }
    }

    // Broadcast grid parameters to all processes
    MPI_Bcast(&NX, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NY, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Lx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Ly, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Create a default configuration and override dt and T.
    SimulationConfig config = default_config();
    config.dt = dt;
    config.T = T;

    // Initialize the simulation with the grid and configuration.
    Simulation *sim = initialize_simulation(grid, config);
    if (!sim) {
        fprintf(stderr, "Error initializing simulation\n");
        if (mpi_data->rank == 0) free_grid(grid);
        free_mpi_data(mpi_data);
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // Distribute initial grid data
    if (mpi_data->rank == 0) {
        distribute_grid(sim, mpi_data);
    }

    int step = 0;
    char vtk_filename[256];
    while (sim->time < sim->T) {
        // Exchange ghost cells before each step
        exchange_ghost_cells(sim, mpi_data);

        // Advance the simulation by one step.
        simulation_step(sim);
        step++;

        // Gather results for output (only on rank 0)
        if (mpi_data->rank == 0 && step % output_interval == 0) {
            gather_grid(sim, mpi_data);
            sprintf(vtk_filename, "output_%04d.vtk", step / output_interval);
            write_vtk(sim, vtk_filename);
        }
    }

    // Cleanup memory.
    free_simulation(sim);
    if (mpi_data->rank == 0) {
        free_grid(grid);
    }
    free_mpi_data(mpi_data);
    
    MPI_Finalize();
    return 0;
}
