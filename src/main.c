// main.c
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "simulation.h"
#include "grid.h"
#include "visualization.h"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Grid & sim params
    const int NX = 100, NY = 100;
    const double Lx = 1.0, Ly = 1.0;
    const double dt = 0.005, T = 5.0;
    const int output_interval = 1;

    // Build grid & config
    Grid *grid = create_grid(NX, NY, Lx, Ly);
    SimulationConfig cfg = default_config();
    cfg.dt = dt;
    cfg.T  = T;

    // Initialize our local stripe
    Simulation *sim = initialize_simulation(grid, cfg, rank, nprocs);

    // Precompute counts/displacements for MPI_Gatherv on rank 0
    int *recvcounts = NULL, *displs = NULL;
    if (rank == 0) {
        recvcounts = malloc(nprocs * sizeof(int));
        displs     = malloc(nprocs * sizeof(int));
        int base = NY / nprocs;
        int rem  = NY % nprocs;
        int offset_rows = 0;
        for (int r = 0; r < nprocs; r++) {
            int rows = base + (r < rem ? 1 : 0);
            recvcounts[r] = rows * NX;
            displs[r]     = offset_rows * NX;
            offset_rows  += rows;
        }
    }

    int step = 0;
    char fname[256];

    while (sim->time < sim->T - 1e-12) {
        // Advance local stripe
        simulation_step(sim);
        step++;

        // Write output every output_interval steps
        if (step % output_interval == 0 && rank == 0) {
            // Gather into a global buffer
            double *global_density = malloc(NX * NY * sizeof(double));
            MPI_Gatherv(
                /* sendbuf: skip ghost‐row 0 */ 
                &sim->density[IX(0,1,NX)],
                /* sendcount */ sim->local_NY * NX,
                MPI_DOUBLE,
                /* recvbuf */ global_density,
                /* recvcounts, displacements */ recvcounts, displs,
                MPI_DOUBLE,
                /* root */ 0,
                MPI_COMM_WORLD
            );

            // Swap in the global buffer for writing
            double *save_ptr = sim->density;
            sim->density = global_density;

            // Write VTK
            sprintf(fname, "output_%04d.vtk", step / output_interval);
            write_vtk(sim, fname);

            // Restore and free
            sim->density = save_ptr;
            free(global_density);
        } else {
            // non-root ranks still must participate in the gather
            MPI_Gatherv(
                &sim->density[IX(0,1,NX)],
                sim->local_NY * NX, MPI_DOUBLE,
                NULL, NULL, NULL, MPI_DOUBLE,
                0, MPI_COMM_WORLD
            );
        }
    }

    // Clean up
    if (rank == 0) {
        free(recvcounts);
        free(displs);
    }
    free_simulation(sim);
    free_grid(grid);
    MPI_Finalize();
    return 0;
}
