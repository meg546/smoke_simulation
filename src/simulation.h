#ifndef SIMULATION_H
#define SIMULATION_H

#include "grid.h"
#include "utils.h"
#include <stdbool.h>
#include <math.h>

// Macro to compute a 1D index for a 2D array stored in row-major order.
#define IX(i, j, NX) ((i) + (NX) * (j))

// -----------------------------------------------------------------------------
// Configuration Structure
// -----------------------------------------------------------------------------

typedef struct {
    double dt;                        // Time step
    double T;                         // Total simulation time
    double density_diffusion;         // Diffusion coefficient for density
    double velocity_viscosity;        // Viscosity for velocity diffusion
    double pressure_tolerance;        // Convergence tolerance for pressure solver
    int    pressure_max_iter;         // Maximum iterations for pressure solver
    double buoyancy_beta;             // Buoyancy coefficient
    double ambient_density;           // Ambient density for buoyancy force
    double wind_u;                    // Horizontal wind velocity
    double wind_v;                    // Vertical wind velocity
    double turbulence_magnitude;      // Turbulence injection magnitude
    int    turbulence_injection_interval; // Steps between turbulence injections
    double vorticity_epsilon;         // Strength of vorticity confinement
    int    smoke_radius;              // Radius for smoke injection
    double smoke_amount;              // Amount of smoke to inject
    int    smoke_pulse_period;
    int    smoke_pulse_duration;
    double velocity_clamp_min;        // Minimum allowed velocity (for clamping)
    double velocity_clamp_max;        // Maximum allowed velocity (for clamping)

    double thermal_diffusivity;
    int ambient_temperature;

} SimulationConfig;

// Returns a default configuration.
SimulationConfig default_config(void);

// -----------------------------------------------------------------------------
// Simulation Structure
// -----------------------------------------------------------------------------

typedef struct Simulation {
    Grid *grid;       // Pointer to the grid structure.
    int NX;           // Number of grid cells in the x-direction.
    int NY;           // Number of grid cells in the y-direction.
    double dt;        // Time step size.
    double T;         // Total simulation time.
    double time;      // Current simulation time.
    
    double *density;  // Smoke density field (size: NX * NY).
    double *u;        // x-component of the velocity field (size: NX * NY).
    double *v;        // y-component of the velocity field (size: NX * NY).
    double *pressure; // Pressure field (size: NX * NY).

    double *temperature;

    SimulationConfig config; // Simulation configuration parameters.
} Simulation;

// -----------------------------------------------------------------------------
// Function Prototypes
// -----------------------------------------------------------------------------

Simulation *initialize_simulation(Grid *grid, SimulationConfig config);
void free_simulation(Simulation *sim);
void advect_density(Simulation *sim);
void advect_temperature(Simulation *sim);
void diffuse_density(Simulation *sim);
void diffuse_velocity(Simulation *sim);
void diffuse_temperature(Simulation *sim);
void solve_pressure(Simulation *sim);
void update_velocity(Simulation *sim);
void apply_buoyancy(Simulation *sim);
void apply_wind(Simulation *sim);
void apply_vorticity_confinement(Simulation *sim);
void inject_turbulence(Simulation *sim);
void inject_smoke(Simulation *sim, int x0, int y0, int radius, double amount);
void inject_heat(Simulation *sim, int x0, int y0, int radius, double heat_amount);
void print_density(const Simulation *sim);
void simulation_step(Simulation *sim);

#endif // SIMULATION_H
