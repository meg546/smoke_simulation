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

    // Create data directory at startup (only rank 0)
    if (rank == 0) {
        system("mkdir -p data");
        // Clear any existing output files
        system("rm -f data/output_*.vtk");
    }
    // Make sure all processes wait until directory is created
    MPI_Barrier(MPI_COMM_WORLD);

    // Grid & sim params
    const int NX = 500, NY = 500;
    const double Lx = 1.0, Ly = 1.0;
    const double dt = 0.0002, T = 1.0;
    const int output_interval = 2;

    // Start total time measurement
    double total_start_time = MPI_Wtime();

    // Build grid & config
    Grid *grid = create_grid(NX, NY, Lx, Ly);
    SimulationConfig cfg = default_config();
    cfg.dt = dt;
    cfg.T  = T;

    // Initialize simulation
    Simulation *sim = initialize_simulation(grid, cfg, rank, nprocs);

    // Barrier to ensure all processes start computation together
    MPI_Barrier(MPI_COMM_WORLD);
    double compute_start_time = MPI_Wtime();

    // Main simulation loop
    int step = 0;
    char fname[256];
    
    while (sim->time < sim->T - 1e-12) {
        simulation_step(sim);
        step++;

        // Write output every output_interval steps
        if (step % output_interval == 0) {
            // Gather density data to rank 0
            double *global_density = NULL;
            if (rank == 0) {
                global_density = malloc(NX * NY * sizeof(double));
            }

            // Gather from all processes to rank 0
            MPI_Gatherv(
                &sim->density[IX(0,1,NX)], 
                sim->local_NY * NX,
                MPI_DOUBLE,
                global_density,
                sim->recvcounts,
                sim->displs,
                MPI_DOUBLE,
                0,
                MPI_COMM_WORLD
            );

            // Write VTK file from rank 0
            if (rank == 0) {
                sprintf(fname, "output_%04d.vtk", step / output_interval);
                write_vtk_global(NX, NY, sim->grid->dx, sim->grid->dy,
                               global_density, fname);
                free(global_density);
            }
        }
    }

    // Ensure all processes finish before stopping timer
    MPI_Barrier(MPI_COMM_WORLD);
    double compute_end_time = MPI_Wtime();
    double total_compute_time = compute_end_time - compute_start_time;

    // Gather timing data from all processes
    double *all_compute_times = NULL;
    if (rank == 0) {
        all_compute_times = (double*)malloc(nprocs * sizeof(double));
    }
    
    MPI_Gather(&total_compute_time, 1, MPI_DOUBLE, 
               all_compute_times, 1, MPI_DOUBLE, 
               0, MPI_COMM_WORLD);

    // Calculate parallel metrics
    if (rank == 0) {
        // Read or write sequential time
        double sequential_time;
        if (nprocs == 1) {
            sequential_time = total_compute_time;
            FILE *fp = fopen("sequential_baseline.txt", "w");
            if (fp) {
                fprintf(fp, "%.6f\n", sequential_time);
                fclose(fp);
            }
        } else {
            FILE *fp = fopen("sequential_baseline.txt", "r");
            if (fp) {
                fscanf(fp, "%lf", &sequential_time);
                fclose(fp);
            } else {
                printf("Warning: Sequential baseline not found. Run with p=1 first.\n");
                sequential_time = -1;
            }
        }

        // Calculate metrics
        sim->metrics.sequential_time = sequential_time;
        sim->metrics.parallel_time = total_compute_time;
        sim->metrics.num_processors = nprocs;

        // Calculate average times per step
        double avg_comp_time = sim->timing.computation_time / sim->timing.step_count;
        double avg_comm_time = sim->timing.communication_time / sim->timing.step_count;
        double avg_pressure_time = sim->timing.pressure_solve_time / sim->timing.step_count;
        double avg_advection_time = sim->timing.advection_time / sim->timing.step_count;
        double avg_total_time = total_compute_time / sim->timing.step_count;

        printf("\n=== Performance Analysis ===\n");
        printf("Number of Processors (p): %d\n", nprocs);
        printf("Problem Size: %d x %d grid\n", sim->NX, sim->NY);
        printf("Total Steps: %d\n", sim->timing.step_count);
        
        printf("\nTotal Times:\n");
        printf("Total Wall Clock Time: %.3f seconds\n", total_compute_time);
        printf("Total Computation Time: %.3f seconds\n", sim->timing.computation_time);
        printf("Total Communication Time: %.3f seconds\n", sim->timing.communication_time);
        
        printf("\nAverage Times per Step:\n");
        printf("Computation: %.3f ms\n", avg_comp_time * 1000.0);
        printf("Communication: %.3f ms\n", avg_comm_time * 1000.0);
        printf("Pressure Solve: %.3f ms\n", avg_pressure_time * 1000.0);
        printf("Advection: %.3f ms\n", avg_advection_time * 1000.0);
        printf("Total Step Time: %.3f ms\n", avg_total_time * 1000.0);

        printf("\nLoad Balance:\n");
        double min_time = all_compute_times[0];
        double max_time = all_compute_times[0];
        double avg_time = all_compute_times[0];
        for (int i = 1; i < nprocs; i++) {
            min_time = fmin(min_time, all_compute_times[i]);
            max_time = fmax(max_time, all_compute_times[i]);
            avg_time += all_compute_times[i];
        }
        avg_time /= nprocs;
        
        printf("Min Process Time: %.3f seconds\n", min_time);
        printf("Max Process Time: %.3f seconds\n", max_time);
        printf("Avg Process Time: %.3f seconds\n", avg_time);
        printf("Load Imbalance: %.2f%%\n",
               ((max_time - min_time) / avg_time) * 100.0);

        if (sequential_time > 0) {
            sim->metrics.speedup = sequential_time / total_compute_time;
            sim->metrics.efficiency = sim->metrics.speedup / nprocs;
            sim->metrics.serial_fraction = (nprocs - sim->metrics.speedup) / 
                                         (sim->metrics.speedup * (nprocs - 1));

            printf("\nParallel Metrics:\n");
            printf("Sequential Time (p=1): %.3f seconds\n", sequential_time);
            printf("Speedup: %.2f\n", sim->metrics.speedup);
            printf("Efficiency: %.2f%%\n", sim->metrics.efficiency * 100.0);
            printf("Estimated Serial Fraction: %.3f\n", sim->metrics.serial_fraction);
            
            printf("\nScaling Analysis:\n");
            printf("Linear Speedup (theoretical): %.2f\n", (double)nprocs);
            printf("Scaling Efficiency: %.2f%%\n",
                   (sim->metrics.speedup / nprocs) * 100.0);

            // --- append CSV with all metrics for this run ---
            const char *outfname = "scaling_results.csv";
            FILE *fp = fopen(outfname, "a");
            if (fp) {
                // if the file was empty, write header
                fseek(fp, 0, SEEK_END);
                if (ftell(fp) == 0) {
                    fprintf(fp,
                        "p,wall_time_s,avg_comp_ms,avg_comm_ms,speedup,efficiency_pct,serial_frac\n");
                }
                // compute averages in ms
                double avg_comp_ms = (sim->timing.computation_time / sim->timing.step_count) * 1e3;
                double avg_comm_ms = (sim->timing.communication_time  / sim->timing.step_count) * 1e3;

                fprintf(fp,
                    "%d,%.6f,%.3f,%.3f,%.3f,%.1f,%.6f\n",
                    nprocs,
                    total_compute_time,
                    avg_comp_ms,
                    avg_comm_ms,
                    sim->metrics.speedup,
                    sim->metrics.efficiency * 100.0,
                    sim->metrics.serial_fraction
                );
                fclose(fp);
            }
        }
    }

    // Cleanup
    if (rank == 0) {
        free(all_compute_times);
    }
    free_simulation(sim);
    free_grid(grid);
    MPI_Finalize();
    return 0;
}
