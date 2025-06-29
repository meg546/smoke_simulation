/*
 * Example Configuration Parameters for Smoke Simulation
 * 
 * This file shows how to modify simulation parameters in main.c
 * Copy these settings to main.c and recompile to use them.
 */

// =============================================================================
// GRID CONFIGURATION
// =============================================================================

// High resolution (slow but detailed)
const int NX = 512, NY = 512;
const double Lx = 1.0, Ly = 1.0;

// Medium resolution (balanced)
// const int NX = 256, NY = 256;
// const double Lx = 1.0, Ly = 1.0;

// Low resolution (fast testing)
// const int NX = 128, NY = 128;
// const double Lx = 1.0, Ly = 1.0;

// =============================================================================
// TIME STEPPING
// =============================================================================

// Short simulation (quick test)
const double dt = 0.003, T = 1.0;

// Medium simulation (good for demos)
// const double dt = 0.003, T = 3.0;

// Long simulation (detailed evolution)
// const double dt = 0.002, T = 10.0;

// =============================================================================
// PHYSICS PARAMETERS (modify in simulation.c default_config())
// =============================================================================

/*
// Weak smoke (subtle effects)
c.density_diffusion = 0.0002;
c.velocity_viscosity = 0.001;
c.buoyancy_beta = 0.01;
c.smoke_amount = 0.2;

// Standard smoke (balanced)
c.density_diffusion = 0.0001;
c.velocity_viscosity = 0.0005;
c.buoyancy_beta = 0.03;
c.smoke_amount = 0.4;

// Strong smoke (dramatic effects)
c.density_diffusion = 0.00005;
c.velocity_viscosity = 0.0002;
c.buoyancy_beta = 0.05;
c.smoke_amount = 0.6;

// High turbulence
c.turbulence_magnitude = 0.05;
c.vorticity_epsilon = 0.30;

// Low turbulence
c.turbulence_magnitude = 0.01;
c.vorticity_epsilon = 0.10;

// Strong wind
c.wind_u = 0.10;
c.wind_v = 0.02;

// No wind
c.wind_u = 0.0;
c.wind_v = 0.0;
*/

// =============================================================================
// OUTPUT CONFIGURATION
// =============================================================================

// High frequency output (every step)
const int output_interval = 1;

// Medium frequency output (every 5 steps)
// const int output_interval = 5;

// Low frequency output (every 10 steps)
// const int output_interval = 10;

// =============================================================================
// USAGE INSTRUCTIONS
// =============================================================================

/*
1. Copy desired parameters to src/main.c
2. Recompile: make clean && make
3. Run simulation: mpirun -n 4 ./smoke_sim
4. Visualize in ParaView using data/output_*.vtk files

For parameter tuning:
- Increase resolution (NX, NY) for more detail
- Decrease dt for numerical stability
- Increase T for longer simulations
- Adjust physics parameters for different effects
*/
