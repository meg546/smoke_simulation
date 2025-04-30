// simulation.h
#ifndef SIMULATION_H
#define SIMULATION_H

#include "grid.h"
#include "utils.h"
#include <mpi.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

// Macro to compute a 1D index for a 2D array stored in row-major order.
#define IX(i,j,stride) ((i) + (stride) * (j))

typedef struct {
    double dt;
    double T;
    double density_diffusion;
    double velocity_viscosity;
    double pressure_tolerance;
    int    pressure_max_iter;
    double buoyancy_beta;
    double ambient_density;
    double wind_u;
    double wind_v;
    double turbulence_magnitude;
    int    turbulence_injection_interval;
    double vorticity_epsilon;
    int    smoke_radius;
    double smoke_amount;
    int    smoke_pulse_period;
    int    smoke_pulse_duration;
    double velocity_clamp_min;
    double velocity_clamp_max;
    double thermal_diffusivity;
    int    ambient_temperature;
} SimulationConfig;

typedef struct SimulationTiming {
    double total_time;
    double computation_time;
    double communication_time;
    double pressure_solve_time;
    double advection_time;
    double diffusion_time;
    double force_time;
    double io_time;
    int step_count;
} SimulationTiming;

typedef struct ParallelMetrics {
    double sequential_time;    // Time for p=1 (baseline)
    double parallel_time;      // Time for current p
    int num_processors;        // Current p value
    double speedup;           // Sequential_time / Parallel_time
    double efficiency;        // Speedup / p
    double serial_fraction;   // Calculated using Amdahl's Law
} ParallelMetrics;

typedef struct Simulation {
    Grid *grid;
    int NX, NY;            // GLOBAL dimensions
    double dt, T, time;
    SimulationConfig config;

    // MPI domain decomposition
    int  rank, nprocs;
    int  local_NY;         // number of interior rows this rank owns
    int  j0;               // global j-index of our first interior row

    // Each array has size NX * (local_NY + 2) to store two ghost‐rows
    double *density, *u, *v, *pressure, *temperature;
    SimulationTiming timing;
    ParallelMetrics metrics;
    int *recvcounts;  // Array of receive counts for MPI_Gatherv
    int *displs;      // Array of displacements for MPI_Gatherv
} Simulation;

SimulationConfig default_config(void);
Simulation *initialize_simulation(Grid *grid, SimulationConfig config,
                                  int rank, int nprocs);
void free_simulation(Simulation *sim);

// core physics on each rank's local stripe
void simulation_step(Simulation *sim);

// I/O (only rank 0 writes)
void write_vtk(Simulation *sim, const char *filename);

#endif // SIMULATION_H
