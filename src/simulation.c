#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include "simulation.h"
#include "grid.h"
#include "utils.h"
#include "mpi_utils.h"

// -----------------------------------------------------------------------------
// Default Configuration
// -----------------------------------------------------------------------------

SimulationConfig default_config() {
    SimulationConfig config;
    config.density_diffusion = 0.0001;
    config.velocity_viscosity = 0.005;
    config.pressure_tolerance = 1e-6;
    config.pressure_max_iter = 1000;
    config.buoyancy_beta = 0.08;
    config.ambient_temperature = 300.0;
    config.ambient_density = 0.01;
    config.wind_u = 4.5;
    config.wind_v = 3.5;
    config.turbulence_magnitude = 0.02;
    config.turbulence_injection_interval = 2;
    config.vorticity_epsilon = 0.2;
    config.smoke_radius = 8;
    config.smoke_amount = 0.4;
    config.smoke_pulse_period = 10;    // Inject smoke every 20 steps.
    config.smoke_pulse_duration = 5;     // But only for 5 steps of each period.
    config.velocity_clamp_min = -5.0;
    config.velocity_clamp_max = 5.0;

    config.thermal_diffusivity = 0.0001;

    return config;
}

// -----------------------------------------------------------------------------
// Memory Management and Initialization
// -----------------------------------------------------------------------------

Simulation *initialize_simulation(Grid *grid, SimulationConfig config)
{
    Simulation *sim = malloc(sizeof(Simulation));
    if (!sim) {
        fprintf(stderr, "Error allocating memory for simulation\n");
        return NULL;
    }

    sim->grid = grid;
    sim->NX = grid->NX;
    sim->NY = grid->NY;
    sim->dt = config.dt;
    sim->T = config.T;
    sim->time = 0.0;
    sim->config = config;  // Save configuration

    int size = sim->NX * sim->NY;

    sim->density = calloc(size, sizeof(double));
    sim->u = calloc(size, sizeof(double));
    sim->v = calloc(size, sizeof(double));
    sim->pressure = calloc(size, sizeof(double));
    sim->temperature = calloc(size, sizeof(double));

    if (!sim->density || !sim->u || !sim->v || !sim->pressure) {
        fprintf(stderr, "Error allocating simulation fields\n");
        free_simulation(sim);
        return NULL;
    }

    // Set the temperature field to the ambient temperature.
    for (int i = 0; i < size; i++) {
        sim->temperature[i] = config.ambient_temperature;
    }

    return sim;
}

void free_simulation(Simulation *sim)
{
    if (!sim) return;
    free(sim->density);
    free(sim->u);
    free(sim->v);
    free(sim->pressure);
    free(sim->temperature);
    free(sim);
}

// -----------------------------------------------------------------------------
// Advection
// -----------------------------------------------------------------------------

void advect_density(Simulation *sim, MPI_Data *mpi_data) {
    int NX = sim->NX;
    int local_ny = mpi_data->local_ny;
    double dt = sim->dt;
    double dx = sim->grid->dx, dy = sim->grid->dy;
    int size = NX * local_ny;

    double *newDensity = malloc(size * sizeof(double));
    if (!newDensity) {
        fprintf(stderr, "Error allocating memory for new density field\n");
        return;
    }

    for (int j = 0; j < local_ny; j++) {
        for (int i = 0; i < NX; i++) {
            int idx = IX(i, j, NX);
            double x = i * dx;
            double y = (j + mpi_data->start_y) * dy;
            double u = sim->u[idx];
            double v = sim->v[idx];

            // Backtrace to find source position.
            double x_src = x - u * dt;
            double y_src = y - v * dt;

            // Clamp source position to grid bounds.
            if (x_src < 0.0) x_src = 0.0;
            else if (x_src > (NX - 1) * dx) x_src = (NX - 1) * dx;
            if (y_src < 0.0) y_src = 0.0;
            else if (y_src > (sim->NY - 1) * dy) y_src = (sim->NY - 1) * dy;

            // Bilinear interpolation.
            int i0 = (int)floor(x_src / dx);
            int j0 = (int)floor(y_src / dy);
            int i1 = (i0 + 1 >= NX) ? NX - 1 : i0 + 1;
            int j1 = (j0 + 1 >= sim->NY) ? sim->NY - 1 : j0 + 1;
            double sx = (x_src / dx) - i0;
            double sy = (y_src / dy) - j0;

            // Handle interpolation across process boundaries
            double d00, d01, d10, d11;
            
            // Get values from local or ghost cells
            if (j0 >= mpi_data->start_y && j0 < mpi_data->start_y + local_ny) {
                d00 = sim->density[IX(i0, j0 - mpi_data->start_y, NX)];
            } else {
                d00 = 0.0; // Use ghost cell value
            }
            
            if (j1 >= mpi_data->start_y && j1 < mpi_data->start_y + local_ny) {
                d01 = sim->density[IX(i0, j1 - mpi_data->start_y, NX)];
            } else {
                d01 = 0.0; // Use ghost cell value
            }
            
            if (j0 >= mpi_data->start_y && j0 < mpi_data->start_y + local_ny) {
                d10 = sim->density[IX(i1, j0 - mpi_data->start_y, NX)];
            } else {
                d10 = 0.0; // Use ghost cell value
            }
            
            if (j1 >= mpi_data->start_y && j1 < mpi_data->start_y + local_ny) {
                d11 = sim->density[IX(i1, j1 - mpi_data->start_y, NX)];
            } else {
                d11 = 0.0; // Use ghost cell value
            }

            double density_interp =
                (1 - sx) * (1 - sy) * d00 +
                sx * (1 - sy) * d10 +
                (1 - sx) * sy * d01 +
                sx * sy * d11;

            newDensity[idx] = density_interp;
        }
    }

    memcpy(sim->density, newDensity, size * sizeof(double));
    free(newDensity);
}

void advect_temperature(Simulation *sim) {
    int NX = sim->NX, NY = sim->NY;
    double dt = sim->dt;
    double dx = sim->grid->dx, dy = sim->grid->dy;
    int size = NX * NY;
    
    double *newTemp = malloc(size * sizeof(double));
    if (!newTemp) {
        fprintf(stderr, "Error allocating memory for temperature advection\n");
        return;
    }
    
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            int idx = IX(i, j, NX);
            double x = i * dx;
            double y = j * dy;
            double u = sim->u[idx];
            double v = sim->v[idx];
            
            // Backtrace the temperature field.
            double x_src = x - u * dt;
            double y_src = y - v * dt;
            
            if (x_src < 0.0) x_src = 0.0;
            else if (x_src > (NX - 1) * dx) x_src = (NX - 1) * dx;
            if (y_src < 0.0) y_src = 0.0;
            else if (y_src > (NY - 1) * dy) y_src = (NY - 1) * dy;
            
            int i0 = (int)floor(x_src / dx);
            int j0 = (int)floor(y_src / dy);
            int i1 = (i0 + 1 >= NX) ? NX - 1 : i0 + 1;
            int j1 = (j0 + 1 >= NY) ? NY - 1 : j0 + 1;
            double sx = (x_src / dx) - i0;
            double sy = (y_src / dy) - j0;
            
            double temp_interp =
                (1 - sx) * (1 - sy) * sim->temperature[IX(i0, j0, NX)] +
                sx * (1 - sy) * sim->temperature[IX(i1, j0, NX)] +
                (1 - sx) * sy * sim->temperature[IX(i0, j1, NX)] +
                sx * sy * sim->temperature[IX(i1, j1, NX)];
            
            newTemp[idx] = temp_interp;
        }
    }
    
    memcpy(sim->temperature, newTemp, size * sizeof(double));
    free(newTemp);
}


// -----------------------------------------------------------------------------
// Diffusion
// -----------------------------------------------------------------------------

void diffuse_density(Simulation *sim)
{
    int NX = sim->NX, NY = sim->NY;
    double dt = sim->dt, dx = sim->grid->dx, dy = sim->grid->dy;
    int size = NX * NY;
    double diff = sim->config.density_diffusion;

    double *newDensity = malloc(size * sizeof(double));
    if (!newDensity) {
        fprintf(stderr, "Error allocating memory for density diffusion\n");
        return;
    }

    for (int j = 1; j < NY - 1; j++) {
        for (int i = 1; i < NX - 1; i++) {
            int idx = IX(i, j, NX);
            double d_center = sim->density[idx];
            double d_left   = sim->density[IX(i - 1, j, NX)];
            double d_right  = sim->density[IX(i + 1, j, NX)];
            double d_down   = sim->density[IX(i, j - 1, NX)];
            double d_up     = sim->density[IX(i, j + 1, NX)];

            double laplacian = (d_left - 2 * d_center + d_right) / (dx * dy)
                             + (d_down - 2 * d_center + d_up) / (dy * dy);
            newDensity[idx] = d_center + dt * diff * laplacian;
        }
    }

    // Copy boundaries (no diffusion applied).
    for (int i = 0; i < NX; i++) {
        newDensity[IX(i, 0, NX)] = sim->density[IX(i, 0, NX)];
        newDensity[IX(i, NY - 1, NX)] = sim->density[IX(i, NY - 1, NX)];
    }
    for (int j = 0; j < NY; j++) {
        newDensity[IX(0, j, NX)] = sim->density[IX(0, j, NX)];
        newDensity[IX(NX - 1, j, NX)] = sim->density[IX(NX - 1, j, NX)];
    }

    memcpy(sim->density, newDensity, size * sizeof(double));
    free(newDensity);
}

void diffuse_velocity(Simulation *sim)
{
    int NX = sim->NX, NY = sim->NY;
    double dt = sim->dt, dx = sim->grid->dx, dy = sim->grid->dy;
    int size = NX * NY;
    double visc = sim->config.velocity_viscosity;

    double *newU = malloc(size * sizeof(double));
    double *newV = malloc(size * sizeof(double));
    if (!newU || !newV) {
        fprintf(stderr, "Error allocating memory for velocity diffusion\n");
        if (newU) free(newU);
        if (newV) free(newV);
        return;
    }

    for (int j = 1; j < NY - 1; j++) {
        for (int i = 1; i < NX - 1; i++) {
            int idx = IX(i, j, NX);
            double u_center = sim->u[idx];
            double laplacian_u = (sim->u[IX(i - 1, j, NX)] - 2 * u_center + sim->u[IX(i + 1, j, NX)]) / (dx * dx)
                               + (sim->u[IX(i, j - 1, NX)] - 2 * u_center + sim->u[IX(i, j + 1, NX)]) / (dy * dy);
            newU[idx] = u_center + dt * visc * laplacian_u;

            double v_center = sim->v[idx];
            double laplacian_v = (sim->v[IX(i - 1, j, NX)] - 2 * v_center + sim->v[IX(i + 1, j, NX)]) / (dx * dx)
                               + (sim->v[IX(i, j - 1, NX)] - 2 * v_center + sim->v[IX(i, j + 1, NX)]) / (dy * dy);
            newV[idx] = v_center + dt * visc * laplacian_v;
        }
    }

    // Copy boundary cells.
    for (int j = 0; j < NY; j++) {
        newU[IX(0, j, NX)] = sim->u[IX(0, j, NX)];
        newU[IX(NX - 1, j, NX)] = sim->u[IX(NX - 1, j, NX)];
        newV[IX(0, j, NX)] = sim->v[IX(0, j, NX)];
        newV[IX(NX - 1, j, NX)] = sim->v[IX(NX - 1, j, NX)];
    }
    for (int i = 0; i < NX; i++) {
        newU[IX(i, 0, NX)] = sim->u[IX(i, 0, NX)];
        newU[IX(i, NY - 1, NX)] = sim->u[IX(i, NY - 1, NX)];
        newV[IX(i, 0, NX)] = sim->v[IX(i, 0, NX)];
        newV[IX(i, NY - 1, NX)] = sim->v[IX(i, NY - 1, NX)];
    }

    memcpy(sim->u, newU, size * sizeof(double));
    memcpy(sim->v, newV, size * sizeof(double));

    free(newU);
    free(newV);
}

void diffuse_temperature(Simulation *sim) {
    int NX = sim->NX, NY = sim->NY;
    double dt = sim->dt, dx = sim->grid->dx, dy = sim->grid->dy;
    int size = NX * NY;
    double alpha = sim->config.thermal_diffusivity;
    
    double *newTemp = malloc(size * sizeof(double));
    if (!newTemp) {
        fprintf(stderr, "Error allocating memory for temperature diffusion\n");
        return;
    }
    
    // Apply a simple finite-difference diffusion scheme.
    for (int j = 1; j < NY - 1; j++) {
        for (int i = 1; i < NX - 1; i++) {
            int idx = IX(i, j, NX);
            double T_center = sim->temperature[idx];
            double T_left = sim->temperature[IX(i - 1, j, NX)];
            double T_right = sim->temperature[IX(i + 1, j, NX)];
            double T_down = sim->temperature[IX(i, j - 1, NX)];
            double T_up = sim->temperature[IX(i, j + 1, NX)];
            
            double laplacian = (T_left - 2 * T_center + T_right) / (dx * dx)
                             + (T_down - 2 * T_center + T_up) / (dy * dy);
            
            newTemp[idx] = T_center + dt * alpha * laplacian;
        }
    }
    
    // Copy boundaries (for simplicity, copying existing temperatures).
    for (int i = 0; i < NX; i++) {
        newTemp[IX(i, 0, NX)] = sim->temperature[IX(i, 0, NX)];
        newTemp[IX(i, NY - 1, NX)] = sim->temperature[IX(i, NY - 1, NX)];
    }
    for (int j = 0; j < NY; j++) {
        newTemp[IX(0, j, NX)] = sim->temperature[IX(0, j, NX)];
        newTemp[IX(NX - 1, j, NX)] = sim->temperature[IX(NX - 1, j, NX)];
    }
    
    memcpy(sim->temperature, newTemp, size * sizeof(double));
    free(newTemp);
}


// -----------------------------------------------------------------------------
// Pressure Solve
// -----------------------------------------------------------------------------

void solve_pressure(Simulation *sim)
{
    int NX = sim->NX, NY = sim->NY;
    double dx = sim->grid->dx, dy = sim->grid->dy, dt = sim->dt;
    int size = NX * NY;
    int max_iter = sim->config.pressure_max_iter;
    double tol = sim->config.pressure_tolerance;

    double *p_new = calloc(size, sizeof(double));
    if (!p_new) {
        fprintf(stderr, "Error allocating memory for pressure solver\n");
        return;
    }

    for (int iter = 0; iter < max_iter; iter++) {
        double max_change = 0.0;
        for (int j = 1; j < NY - 1; j++) {
            for (int i = 1; i < NX - 1; i++) {
                int idx = IX(i, j, NX);
                double divergence = ((sim->u[IX(i+1, j, NX)] - sim->u[IX(i-1, j, NX)]) / (2 * dx))
                                  + ((sim->v[IX(i, j+1, NX)] - sim->v[IX(i, j-1, NX)]) / (2 * dy));
                double new_val = (sim->pressure[IX(i-1, j, NX)] +
                                  sim->pressure[IX(i+1, j, NX)] +
                                  sim->pressure[IX(i, j-1, NX)] +
                                  sim->pressure[IX(i, j+1, NX)] -
                                  (dx * dx / dt) * divergence) / 4.0;
                max_change = fmax(max_change, fabs(new_val - sim->pressure[idx]));
                p_new[idx] = new_val;
            }
        }
        for (int j = 1; j < NY - 1; j++)
            for (int i = 1; i < NX - 1; i++)
                sim->pressure[IX(i, j, NX)] = p_new[IX(i, j, NX)];

        // Enforce Neumann boundary conditions.
        for (int j = 0; j < NY; j++) {
            sim->pressure[IX(0, j, NX)] = sim->pressure[IX(1, j, NX)];
            sim->pressure[IX(NX-1, j, NX)] = sim->pressure[IX(NX-2, j, NX)];
        }
        for (int i = 0; i < NX; i++) {
            sim->pressure[IX(i, 0, NX)] = sim->pressure[IX(i, 1, NX)];
            sim->pressure[IX(i, NY-1, NX)] = sim->pressure[IX(i, NY-2, NX)];
        }
        if (max_change < tol)
            break;
    }
    free(p_new);
}

// -----------------------------------------------------------------------------
// Velocity Update and Clamping
// -----------------------------------------------------------------------------

void update_velocity(Simulation *sim)
{
    int NX = sim->NX, NY = sim->NY;
    double dt = sim->dt, dx = sim->grid->dx, dy = sim->grid->dy;
    int size = NX * NY;

    double *newU = malloc(size * sizeof(double));
    double *newV = malloc(size * sizeof(double));
    if (!newU || !newV) {
        fprintf(stderr, "Error allocating memory in update_velocity\n");
        if (newU) free(newU);
        if (newV) free(newV);
        return;
    }

    for (int j = 1; j < NY - 1; j++) {
        for (int i = 1; i < NX - 1; i++) {
            int idx = IX(i, j, NX);
            double dpdx = (sim->pressure[IX(i+1, j, NX)] - sim->pressure[IX(i-1, j, NX)]) / (2.0 * dx);
            double dpdy = (sim->pressure[IX(i, j+1, NX)] - sim->pressure[IX(i, j-1, NX)]) / (2.0 * dy);
            double u_new = sim->u[idx] - dt * dpdx;
            double v_new = sim->v[idx] - dt * dpdy;
            // Clamp the updated velocities.
            u_new = clamp(u_new, sim->config.velocity_clamp_min, sim->config.velocity_clamp_max);
            v_new = clamp(v_new, sim->config.velocity_clamp_min, sim->config.velocity_clamp_max);
            newU[idx] = u_new;
            newV[idx] = v_new;
        }
    }

    // Enforce no-slip boundary conditions.
    for (int j = 0; j < NY; j++) {
        newU[IX(0, j, NX)] = 0.0;
        newU[IX(NX-1, j, NX)] = 0.0;
        newV[IX(0, j, NX)] = 0.0;
        newV[IX(NX-1, j, NX)] = 0.0;
    }
    for (int i = 0; i < NX; i++) {
        newU[IX(i, 0, NX)] = 0.0;
        newU[IX(i, NY-1, NX)] = 0.0;
        newV[IX(i, 0, NX)] = 0.0;
        newV[IX(i, NY-1, NX)] = 0.0;
    }

    memcpy(sim->u, newU, size * sizeof(double));
    memcpy(sim->v, newV, size * sizeof(double));

    free(newU);
    free(newV);
}

// -----------------------------------------------------------------------------
// External Forces
// -----------------------------------------------------------------------------

void apply_buoyancy(Simulation *sim) {
    int NX = sim->NX, NY = sim->NY;
    double dt = sim->dt;
    double ambientTemp = sim->config.ambient_temperature;
    double beta = sim->config.buoyancy_beta; 

    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            int idx = IX(i, j, NX);
            double buoyancy = beta * (sim->temperature[idx] - ambientTemp);
            sim->v[idx] += dt * buoyancy;
        }
    }
}


void apply_wind(Simulation *sim) {
    int NX = sim->NX, NY = sim->NY;
    double dt = sim->dt;
    double wind_u = sim->config.wind_u, wind_v = sim->config.wind_v;

    for (int j = 1; j < NY - 1; j++) {
        for (int i = 1; i < NX - 1; i++) {
            int idx = IX(i, j, NX);
            sim->u[idx] += wind_u * dt;
            sim->v[idx] += wind_v * dt;
        }
    }
}

void apply_vorticity_confinement(Simulation *sim)
{
    int NX = sim->NX, NY = sim->NY;
    double dx = sim->grid->dx, dy = sim->grid->dy;
    double epsilon = sim->config.vorticity_epsilon;

    double *curl = calloc(NX * NY, sizeof(double));
    if (!curl) return;

    for (int j = 1; j < NY - 1; j++) {
        for (int i = 1; i < NX - 1; i++) {
            int idx = IX(i, j, NX);
            double dv_dx = (sim->v[IX(i+1, j, NX)] - sim->v[IX(i-1, j, NX)]) / (2 * dx);
            double du_dy = (sim->u[IX(i, j+1, NX)] - sim->u[IX(i, j-1, NX)]) / (2 * dy);
            curl[idx] = dv_dx - du_dy;
        }
    }

    for (int j = 2; j < NY - 2; j++) {
        for (int i = 2; i < NX - 2; i++) {
            int idx = IX(i, j, NX);
            double curl_L = fabs(curl[IX(i-1, j, NX)]);
            double curl_R = fabs(curl[IX(i+1, j, NX)]);
            double curl_B = fabs(curl[IX(i, j-1, NX)]);
            double curl_T = fabs(curl[IX(i, j+1, NX)]);
            double N_x = (curl_R - curl_L) / (2 * dx);
            double N_y = (curl_T - curl_B) / (2 * dy);
            double length = sqrt(N_x * N_x + N_y * N_y) + 1e-5;
            N_x /= length;
            N_y /= length;
            sim->u[idx] += epsilon * (dy * N_y * curl[idx]);
            sim->v[idx] -= epsilon * (dx * N_x * curl[idx]);
        }
    }
    free(curl);
}

void inject_turbulence(Simulation *sim) {
    int NX = sim->NX, NY = sim->NY;
    double magnitude = sim->config.turbulence_magnitude;

    for (int j = 1; j < NY - 1; j++) {
        for (int i = 1; i < NX - 1; i++) {
            int idx = IX(i, j, NX);
            sim->u[idx] += magnitude * ((rand() / (double)RAND_MAX) - 0.5);
            sim->v[idx] += magnitude * ((rand() / (double)RAND_MAX) - 0.5);
        }
    }
}

void inject_smoke(Simulation *sim, int x0, int y0, int radius, double amount)
{
    int NX = sim->NX, NY = sim->NY;
    double sigma = radius / 2.0;  // Standard deviation for Gaussian falloff.

    for (int j = y0 - radius; j <= y0 + radius; j++) {
        for (int i = x0 - radius; i <= x0 + radius; i++) {
            if (i < 0 || i >= NX || j < 0 || j >= NY)
                continue;
            int dx = i - x0, dy = j - y0;
            double dist = sqrt(dx * dx + dy * dy);
            if (dist <= radius) {
                int idx = IX(i, j, NX);
                // Gaussian falloff: maximum at the center decays smoothly.
                double factor = exp(- (dist * dist) / (2 * sigma * sigma));
                sim->density[idx] += amount * factor;
            }
        }
    }
}

void inject_heat(Simulation *sim, int x0, int y0, int radius, double heat_amount)
{
    int NX = sim->NX, NY = sim->NY;
    double sigma = radius / 2.0;  // For a Gaussian falloff in temperature

    for (int j = y0 - radius; j <= y0 + radius; j++) {
        for (int i = x0 - radius; i <= x0 + radius; i++) {
            if (i < 0 || i >= NX || j < 0 || j >= NY)
                continue;
            int dx = i - x0, dy = j - y0;
            double dist = sqrt(dx * dx + dy * dy);
            if (dist <= radius) {
                int idx = IX(i, j, NX);
                double factor = exp(- (dist * dist) / (2 * sigma * sigma));
                sim->temperature[idx] += heat_amount * factor;
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Debugging and Output
// -----------------------------------------------------------------------------

void print_density(const Simulation *sim) {
    int NX = sim->NX, NY = sim->NY;
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            int idx = IX(i, j, NX);
            printf("%.2f ", sim->density[idx]);
        }
        printf("\n");
    }
}

// -----------------------------------------------------------------------------
// Main Simulation Step
// -----------------------------------------------------------------------------

void simulation_step(Simulation *sim, MPI_Data *mpi_data)
{
    // Exchange ghost cells before each operation
    exchange_ghost_cells(sim, mpi_data);
    advect_density(sim, mpi_data);
    
    exchange_ghost_cells(sim, mpi_data);
    diffuse_density(sim);
    
    exchange_ghost_cells(sim, mpi_data);
    advect_velocity(sim, mpi_data);
    
    exchange_ghost_cells(sim, mpi_data);
    diffuse_velocity(sim);
    
    exchange_ghost_cells(sim, mpi_data);
    solve_pressure(sim);
    
    exchange_ghost_cells(sim, mpi_data);
    update_velocity(sim);
    
    exchange_ghost_cells(sim, mpi_data);
    apply_buoyancy(sim);
    
    exchange_ghost_cells(sim, mpi_data);
    apply_wind(sim);
    
    exchange_ghost_cells(sim, mpi_data);
    apply_vorticity_confinement(sim);
    
    sim->time += sim->dt;
    printf("Completed simulation step. New time = %f\n", sim->time);
}
