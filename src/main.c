#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simulation.h"
#include "grid.h"
#include "visualization.h"

int main(void) {
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

    // Create the simulation grid.
    Grid *grid = create_grid(NX, NY, Lx, Ly);
    if (!grid) {
        fprintf(stderr, "Error creating grid\n");
        return EXIT_FAILURE;
    }

    // Create a default configuration and override dt and T.
    SimulationConfig config = default_config();
    config.dt = dt;
    config.T = T;

    // Initialize the simulation with the grid and configuration.
    Simulation *sim = initialize_simulation(grid, config);
    if (!sim) {
        fprintf(stderr, "Error initializing simulation\n");
        free_grid(grid);
        return EXIT_FAILURE;
    }

    int step = 0;
    char vtk_filename[256];
    while (sim->time < sim->T) {
        // Only output every 'output_interval' steps.
        if (step % output_interval == 0) {
            sprintf(vtk_filename, "output_%04d.vtk", step / output_interval);
            write_vtk(sim, vtk_filename);
        }
        // Advance the simulation by one step.
        simulation_step(sim);
        step++;
    }

    // Cleanup memory.
    free_simulation(sim);
    free_grid(grid);
    return 0;
}
