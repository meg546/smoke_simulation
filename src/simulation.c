#include "simulation.h"
#include "grid.h"
#include "utils.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// -----------------------------------------------------------------------------
// Exchange our two ghost-rows (j=0 and j=local_NY+1) with up/down neighbors
// -----------------------------------------------------------------------------
static void sim_exchange_ghost_rows(Simulation *sim, double *field) {
    MPI_Status status;
    int NX = sim->NX;
    int M  = sim->local_NY;

    // send our first real row j=1 ↑ to rank-1 → receive into j=0
    if (sim->rank > 0) {
        MPI_Sendrecv(
            &field[IX(0,1,NX)],     NX, MPI_DOUBLE, sim->rank-1, 100,
            &field[IX(0,0,NX)],     NX, MPI_DOUBLE, sim->rank-1, 101,
            MPI_COMM_WORLD, &status
        );
    }
    // send our last real row j=M ↓ to rank+1 → receive into j=M+1
    if (sim->rank < sim->nprocs-1) {
        MPI_Sendrecv(
            &field[IX(0,M,NX)],     NX, MPI_DOUBLE, sim->rank+1, 101,
            &field[IX(0,M+1,NX)],   NX, MPI_DOUBLE, sim->rank+1, 100,
            MPI_COMM_WORLD, &status
        );
    }
}

// -----------------------------------------------------------------------------
// Enforce no‐flux (Neumann) at the global top/bottom for scalars
//   density and temperature ghost rows mirror their adjacent interior rows
// -----------------------------------------------------------------------------
static void enforce_global_scalar_bc(Simulation *sim, double *field) {
    int NX = sim->NX;
    int M  = sim->local_NY;
    // bottom wall
    if (sim->rank == 0) {
        for (int i = 0; i < NX; i++) {
            field[IX(i,0,    NX)] = field[IX(i,1,    NX)];
        }
    }
    // top wall
    if (sim->rank == sim->nprocs - 1) {
        for (int i = 0; i < NX; i++) {
            field[IX(i,M+1,NX)] = field[IX(i,M,   NX)];
        }
    }
}

// -----------------------------------------------------------------------------
// Default simulation parameters
// -----------------------------------------------------------------------------
SimulationConfig default_config(void) {
    SimulationConfig c;
    c.dt                             = 0.01;    // Keep timestep
    c.T                              = 5.0;
    c.density_diffusion              = 0.00015;  // Slightly increased for better spread
    c.velocity_viscosity             = 0.003;    // Reduced to allow more natural flow
    c.pressure_tolerance             = 1e-4;     // Keep this value
    c.pressure_max_iter              = 50;       // Keep this value
    c.buoyancy_beta                  = 0.05;     // Reduced for more gentle rising
    c.ambient_density                = 0.01;
    c.wind_u                         = 0.1;      // Slight horizontal drift
    c.wind_v                         = 0.0;      // No vertical wind
    c.turbulence_magnitude           = 0.015;    // Increased for more natural movement
    c.turbulence_injection_interval  = 2;        // More frequent turbulence
    c.vorticity_epsilon             = 0.15;     // Increased vorticity for better swirling
    c.smoke_radius                  = 6;        // Smaller initial radius
    c.smoke_amount                  = 0.3;      // Reduced density for better detail
    c.smoke_pulse_period           = 8;        // More frequent pulses
    c.smoke_pulse_duration         = 4;        // Shorter pulses
    c.velocity_clamp_min           = -3.0;     // Reduced clamping for smoother motion
    c.velocity_clamp_max           = 3.0;
    c.thermal_diffusivity         = 0.00012;   // Slightly increased heat diffusion
    c.ambient_temperature         = 300.0;
    return c;
}

// -----------------------------------------------------------------------------
// MPI-aware init/cleanup: split NY into stripes with 2 ghost‐rows each
// -----------------------------------------------------------------------------
Simulation *initialize_simulation(Grid *grid,
                                  SimulationConfig cfg,
                                  int rank, int nprocs)
{
    Simulation *sim = malloc(sizeof(*sim));
    sim->grid      = grid;
    sim->NX        = grid->NX;
    sim->NY        = grid->NY;     // total global rows
    sim->dt        = cfg.dt;
    sim->T         = cfg.T;
    sim->time      = 0.0;
    sim->config    = cfg;
    sim->rank      = rank;
    sim->nprocs    = nprocs;

    // split NY into nearly‐equal stripes
    int base = grid->NY / nprocs;
    int rem  = grid->NY % nprocs;
    sim->local_NY = base + (rank < rem ? 1 : 0);
    // global start‐row index for our interior j=1
    sim->j0       = rank * base + (rank < rem ? rank : rem);

    // allocate NX × (local_NY+2) for each field
    int total = sim->NX * (sim->local_NY + 2);
    sim->density     = calloc(total, sizeof(double));
    sim->u           = calloc(total, sizeof(double));
    sim->v           = calloc(total, sizeof(double));
    sim->pressure    = calloc(total, sizeof(double));
    sim->temperature = calloc(total, sizeof(double));

    // initialize interior to ambient temperature
    for (int j = 1; j <= sim->local_NY; ++j) {
        for (int i = 0; i < sim->NX; ++i) {
            sim->temperature[IX(i,j,sim->NX)] = cfg.ambient_temperature;
        }
    }

    // Initialize timing structure
    sim->timing.total_time = 0.0;
    sim->timing.computation_time = 0.0;
    sim->timing.communication_time = 0.0;
    sim->timing.pressure_solve_time = 0.0;
    sim->timing.advection_time = 0.0;
    sim->timing.diffusion_time = 0.0;
    sim->timing.force_time = 0.0;
    sim->timing.io_time = 0.0;
    sim->timing.step_count = 0;

    // Initialize recvcounts and displs for MPI_Gatherv
    if (rank == 0) {
        sim->recvcounts = malloc(nprocs * sizeof(int));
        sim->displs = malloc(nprocs * sizeof(int));
        int base = grid->NY / nprocs;
        int rem = grid->NY % nprocs;
        int offset = 0;
        for (int r = 0; r < nprocs; r++) {
            int rows = base + (r < rem ? 1 : 0);
            sim->recvcounts[r] = rows * grid->NX;
            sim->displs[r] = offset * grid->NX;
            offset += rows;
        }
    } else {
        sim->recvcounts = NULL;
        sim->displs = NULL;
    }

    return sim;
}

void free_simulation(Simulation *sim) {
    if (!sim) return;
    free(sim->density);
    free(sim->u);
    free(sim->v);
    free(sim->pressure);
    free(sim->temperature);
    free(sim->recvcounts);
    free(sim->displs);
    free(sim);
}

// -----------------------------------------------------------------------------
// Advection kernels (semi-Lagrangian, use pre-exchanged ghost rows)
// -----------------------------------------------------------------------------
void advect_density(Simulation *sim) {
  int NX = sim->NX, M = sim->local_NY;
  double dt = sim->dt, dx = sim->grid->dx, dy = sim->grid->dy;
  double *old = sim->density;
  double *newD = malloc(NX*(M+2)*sizeof(double));

  memcpy(newD, old, NX*(M+2)*sizeof(double));

  for (int j = 1; j <= M; ++j) {
      double y = (sim->j0 + j - 1) * dy;
      for (int i = 0; i < NX; ++i) {
          int idx = IX(i,j,NX);
          double x = i*dx;
          // FIX: Use velocity components from velocity fields, not density!
          double u = sim->u[idx], v = sim->v[idx];
          double xs = fmax(0.0, fmin(x - u*dt, (NX-1)*dx));
          double ys = fmax(0.0, fmin(y - v*dt, (sim->NY-1)*dy));

          // Rest of the function remains the same
          int i0  = (int)floor(xs/dx),
              j0g = (int)floor(ys/dy);
          int i1  = (i0+1 < NX ? i0+1 : NX-1),
              j1g = (j0g+1 < sim->NY ? j0g+1 : sim->NY-1);

          int lj0 = j0g - sim->j0 + 1;
          int lj1 = j1g - sim->j0 + 1;
          lj0 = fmax(0, fmin(lj0, M+1));
          lj1 = fmax(0, fmin(lj1, M+1));

          double sx = xs/dx - i0,
                 sy = ys/dy - j0g;

          double d00 = old[IX(i0,lj0,NX)],
                 d10 = old[IX(i1,lj0,NX)],
                 d01 = old[IX(i0,lj1,NX)],
                 d11 = old[IX(i1,lj1,NX)];

          newD[idx] = (1-sx)*(1-sy)*d00
                     + sx*(1-sy)*d10
                     + (1-sx)*sy  *d01
                     + sx*sy      *d11;
      }
  }

  memcpy(sim->density, newD, NX*(M+2)*sizeof(double));
  free(newD);
}

void advect_temperature(Simulation *sim) {
    int NX = sim->NX, M = sim->local_NY;
    double dt = sim->dt, dx = sim->grid->dx, dy = sim->grid->dy;
    double *old = sim->temperature;
    double *newT = malloc(NX*(M+2)*sizeof(double));

    memcpy(newT, old, NX*(M+2)*sizeof(double));

    for (int j = 1; j <= M; ++j) {
        double y = (sim->j0 + j - 1) * dy;
        for (int i = 0; i < NX; ++i) {
            int idx = IX(i,j,NX);
            double x = i*dx;
            double u = sim->u[idx], v = sim->v[idx];
            double xs = fmax(0.0, fmin(x - u*dt, (NX-1)*dx));
            double ys = fmax(0.0, fmin(y - v*dt, (sim->NY-1)*dy));

            int i0  = (int)floor(xs/dx),
                j0g = (int)floor(ys/dy);
            int i1  = (i0+1 < NX ? i0+1 : NX-1),
                j1g = (j0g+1 < sim->NY ? j0g+1 : sim->NY-1);

            int lj0 = j0g - sim->j0 + 1;
            int lj1 = j1g - sim->j0 + 1;
            lj0 = fmax(0, fmin(lj0, M+1));
            lj1 = fmax(0, fmin(lj1, M+1));

            double sx = xs/dx - i0,
                   sy = ys/dy - j0g;

            double t00 = old[IX(i0,lj0,NX)],
                   t10 = old[IX(i1,lj0,NX)],
                   t01 = old[IX(i0,lj1,NX)],
                   t11 = old[IX(i1,lj1,NX)];

            newT[idx] = (1-sx)*(1-sy)*t00
                       + sx*(1-sy)*t10
                       + (1-sx)*sy  *t01
                       + sx*sy      *t11;
        }
    }

    memcpy(sim->temperature, newT, NX*(M+2)*sizeof(double));
    free(newT);
}

// ----------------------------------------------------------------------------
// Density diffusion (ghosts set by caller)
// ----------------------------------------------------------------------------
void diffuse_density(Simulation *sim) {
    int NX = sim->NX, M = sim->local_NY;
    double dt = sim->dt,
           dx = sim->grid->dx,
           dy = sim->grid->dy,
           diff = sim->config.density_diffusion;
    double *newD = malloc(NX*(M+2)*sizeof(double));

    // interior
    for(int j=1;j<=M;++j){
      for(int i=1;i<NX-1;++i){
        int idx = IX(i,j,NX);
        double c = sim->density[idx];
        double lap = (sim->density[IX(i-1,j,NX)] - 2*c + sim->density[IX(i+1,j,NX)])/(dx*dx)
                   + (sim->density[IX(i,j-1,NX)] - 2*c + sim->density[IX(i,j+1,NX)])/(dy*dy);
        newD[idx] = c + dt*diff*lap;
      }
    }
    // walls
    for(int j=1;j<=M;++j){
      newD[IX(0,j,NX)]    = sim->density[IX(0,j,NX)];
      newD[IX(NX-1,j,NX)] = sim->density[IX(NX-1,j,NX)];
    }
    if(sim->rank==0)
      for(int i=0;i<NX;++i) newD[IX(i,1,NX)] = sim->density[IX(i,1,NX)];
    if(sim->rank==sim->nprocs-1)
      for(int i=0;i<NX;++i) newD[IX(i,M,NX)] = sim->density[IX(i,M,NX)];

    memcpy(sim->density, newD, NX*(M+2)*sizeof(double));
    free(newD);
}

// ----------------------------------------------------------------------------
// Velocity diffusion
// ----------------------------------------------------------------------------
void diffuse_velocity(Simulation *sim) {
    int NX = sim->NX, M = sim->local_NY;
    double dt = sim->dt,
           dx = sim->grid->dx,
           dy = sim->grid->dy,
           visc = sim->config.velocity_viscosity;
    double *u2 = malloc(NX*(M+2)*sizeof(double));
    double *v2 = malloc(NX*(M+2)*sizeof(double));

    for(int j=1;j<=M;++j){
      for(int i=1;i<NX-1;++i){
        int idx = IX(i,j,NX);
        double uc = sim->u[idx], vc = sim->v[idx];
        double lapu = (sim->u[IX(i-1,j,NX)] - 2*uc + sim->u[IX(i+1,j,NX)])/(dx*dx)
                    + (sim->u[IX(i,j-1,NX)] - 2*uc + sim->u[IX(i,j+1,NX)])/(dy*dy);
        double lapv = (sim->v[IX(i-1,j,NX)] - 2*vc + sim->v[IX(i+1,j,NX)])/(dx*dx)
                    + (sim->v[IX(i,j-1,NX)] - 2*vc + sim->v[IX(i,j+1,NX)])/(dy*dy);
        u2[idx] = uc + dt*visc*lapu;
        v2[idx] = vc + dt*visc*lapv;
      }
    }
    // no-slip walls
    for(int j=1;j<=M;++j){
      u2[IX(0,j,NX)]    = 0.0;  v2[IX(0,j,NX)]    = 0.0;
      u2[IX(NX-1,j,NX)] = 0.0;  v2[IX(NX-1,j,NX)] = 0.0;
    }
    if(sim->rank==0)
      for(int i=0;i<NX;++i) u2[IX(i,1,NX)] = v2[IX(i,1,NX)] = 0.0;
    if(sim->rank==sim->nprocs-1)
      for(int i=0;i<NX;++i) u2[IX(i,M,NX)] = v2[IX(i,M,NX)] = 0.0;

    memcpy(sim->u, u2, NX*(M+2)*sizeof(double));
    memcpy(sim->v, v2, NX*(M+2)*sizeof(double));
    free(u2); free(v2);
}

// ----------------------------------------------------------------------------
// Temperature diffusion
// ----------------------------------------------------------------------------
void diffuse_temperature(Simulation *sim) {
    int NX = sim->NX, M = sim->local_NY;
    double dt = sim->dt,
           dx = sim->grid->dx,
           dy = sim->grid->dy,
           alpha = sim->config.thermal_diffusivity;
    double *t2 = malloc(NX*(M+2)*sizeof(double));

    for(int j=1;j<=M;++j){
      for(int i=1;i<NX-1;++i){
        int idx = IX(i,j,NX);
        double tc = sim->temperature[idx];
        double lap = (sim->temperature[IX(i-1,j,NX)] - 2*tc + sim->temperature[IX(i+1,j,NX)])/(dx*dx)
                   + (sim->temperature[IX(i,j-1,NX)] - 2*tc + sim->temperature[IX(i,j+1,NX)])/(dy*dy);
        t2[idx] = tc + dt*alpha*lap;
      }
    }
    // insulated walls
    for(int j=1;j<=M;++j){
      t2[IX(0,j,NX)]    = sim->temperature[IX(0,j,NX)];
      t2[IX(NX-1,j,NX)] = sim->temperature[IX(NX-1,j,NX)];
    }
    if(sim->rank==0)
      for(int i=0;i<NX;++i) t2[IX(i,1,NX)] = sim->temperature[IX(i,1,NX)];
    if(sim->rank==sim->nprocs-1)
      for(int i=0;i<NX;++i) t2[IX(i,M,NX)] = sim->temperature[IX(i,M,NX)];

    memcpy(sim->temperature, t2, NX*(M+2)*sizeof(double));
    free(t2);
}

// ----------------------------------------------------------------------------
// Pressure solve (Jacobi) with Neumann walls implied by ghost‐rows
// ----------------------------------------------------------------------------
void solve_pressure(Simulation *sim) {
    int NX = sim->NX, M = sim->local_NY;
    double dx = sim->grid->dx,
           dy = sim->grid->dy,
           dt = sim->dt,
           tol = sim->config.pressure_tolerance;
    int maxit = sim->config.pressure_max_iter;

    double *p2 = calloc(NX*(M+2), sizeof(double));
    
    // Reset pressure to help prevent artifacts
    memset(sim->pressure, 0, NX*(M+2)*sizeof(double));
    
    for(int it=0; it<maxit; ++it){
        double maxchg=0;
        for(int j=1;j<=M;++j){
            for(int i=1;i<NX-1;++i){
                int idx = IX(i,j,NX);
                
                // Calculate divergence using adjacent cells
                double div = (
                    (sim->u[IX(i+1,j,NX)] - sim->u[IX(i-1,j,NX)])/(2.0*dx) +
                    (sim->v[IX(i,j+1,NX)] - sim->v[IX(i,j-1,NX)])/(2.0*dy)
                );
                
                // Standard 5-point Laplacian stencil
                double p_left   = sim->pressure[IX(i-1,j,NX)];
                double p_right  = sim->pressure[IX(i+1,j,NX)];
                double p_bottom = sim->pressure[IX(i,j-1,NX)];
                double p_top    = sim->pressure[IX(i,j+1,NX)];
                double p_center = sim->pressure[idx];
                
                // Gauss-Seidel iteration
                double new_p = (
                    (p_left + p_right)/(dx*dx) +
                    (p_bottom + p_top)/(dy*dy) -
                    div/dt
                ) / (2.0/(dx*dx) + 2.0/(dy*dy));
                
                maxchg = fmax(maxchg, fabs(new_p - p_center));
                p2[idx] = new_p;
            }
        }
        
        // Update pressure field
        memcpy(&sim->pressure[IX(1,1,NX)],
               &p2[IX(1,1,NX)],
               (NX-2)*M*sizeof(double));
        
        // Apply boundary conditions to pressure
        for(int j=1;j<=M;++j){
            sim->pressure[IX(0,j,NX)] = sim->pressure[IX(1,j,NX)];
            sim->pressure[IX(NX-1,j,NX)] = sim->pressure[IX(NX-2,j,NX)];
        }
        
        if(maxchg < tol) break;
    }
    free(p2);
}

// ----------------------------------------------------------------------------
// Velocity update from pressure gradient
// ----------------------------------------------------------------------------
void update_velocity(Simulation *sim) {
    int NX = sim->NX, M = sim->local_NY;
    double dx = sim->grid->dx,
           dy = sim->grid->dy,
           dt = sim->dt,
           vmin = sim->config.velocity_clamp_min,
           vmax = sim->config.velocity_clamp_max;
    double *u2 = malloc(NX*(M+2)*sizeof(double));
    double *v2 = malloc(NX*(M+2)*sizeof(double));

    // Update velocities using central differences for pressure
    for(int j=1;j<=M;++j){
        for(int i=1;i<NX-1;++i){
            int idx = IX(i,j,NX);
            
            // Compute pressure gradients
            double dpdx = (sim->pressure[IX(i+1,j,NX)] - sim->pressure[IX(i-1,j,NX)])/(2.0*dx);
            double dpdy = (sim->pressure[IX(i,j+1,NX)] - sim->pressure[IX(i,j-1,NX)])/(2.0*dy);
            
            // Update velocities
            u2[idx] = sim->u[idx];
            v2[idx] = sim->v[idx];
            
            // Apply pressure correction
            u2[idx] -= dt * dpdx;
            v2[idx] -= dt * dpdy;
            
            // Apply velocity limits
            u2[idx] = clamp(u2[idx], vmin, vmax);
            v2[idx] = clamp(v2[idx], vmin, vmax);
        }
    }

    // Enforce boundary conditions
    for(int j=1;j<=M;++j){
        u2[IX(0,j,NX)] = 0.0;
        v2[IX(0,j,NX)] = 0.0;
        u2[IX(NX-1,j,NX)] = 0.0;
        v2[IX(NX-1,j,NX)] = 0.0;
    }
    
    if(sim->rank==0)
        for(int i=0;i<NX;++i) {
            u2[IX(i,1,NX)] = 0.0;
            v2[IX(i,1,NX)] = 0.0;
        }
    if(sim->rank==sim->nprocs-1)
        for(int i=0;i<NX;++i) {
            u2[IX(i,M,NX)] = 0.0;
            v2[IX(i,M,NX)] = 0.0;
        }

    memcpy(sim->u, u2, NX*(M+2)*sizeof(double));
    memcpy(sim->v, v2, NX*(M+2)*sizeof(double));
    free(u2); free(v2);
}

// ----------------------------------------------------------------------------
// External forces & injection
// ----------------------------------------------------------------------------
void apply_buoyancy(Simulation *sim) {
    int NX = sim->NX, M = sim->local_NY;
    double dt = sim->dt,
           beta = sim->config.buoyancy_beta,
           Tamb = sim->config.ambient_temperature;
    for(int j=1;j<=M;++j)
      for(int i=0;i<NX;++i)
        sim->v[IX(i,j,NX)] += dt * beta * (sim->temperature[IX(i,j,NX)] - Tamb);
}

void apply_wind(Simulation *sim) {
    int NX = sim->NX, M = sim->local_NY;
    double dt = sim->dt,
           wu = sim->config.wind_u,
           wv = sim->config.wind_v;
    for(int j=1;j<=M;++j)
      for(int i=1;i<NX-1;++i){
        sim->u[IX(i,j,NX)] += dt * wu;
        sim->v[IX(i,j,NX)] += dt * wv;
      }
}

void apply_vorticity_confinement(Simulation *sim) {
    int NX = sim->NX, M = sim->local_NY;
    double dx = sim->grid->dx,
           dy = sim->grid->dy,
           eps= sim->config.vorticity_epsilon;
    double *curl = calloc(NX*(M+2), sizeof(double));

    // compute curl
    for(int j=1;j<=M;++j)
      for(int i=1;i<NX-1;++i){
        int idx = IX(i,j,NX);
        double dvdx = (sim->v[IX(i+1,j,NX)] - sim->v[IX(i-1,j,NX)])/(2*dx);
        double dudy = (sim->u[IX(i,j+1,NX)] - sim->u[IX(i,j-1,NX)])/(2*dy);
        curl[idx] = dvdx - dudy;
      }

    // apply confinement
    for(int j=2;j<M;++j)
      for(int i=2;i<NX-2;++i){
        int idx = IX(i,j,NX);
        double L=fabs(curl[IX(i-1,j,NX)]), R=fabs(curl[IX(i+1,j,NX)]);
        double B=fabs(curl[IX(i,j-1,NX)]), T=fabs(curl[IX(i,j+1,NX)]);
        double Nx_ = (R-L)/(2*dx), Ny_ = (T-B)/(2*dy);
        double norm = sqrt(Nx_*Nx_+Ny_*Ny_) + 1e-5;
        Nx_/=norm; Ny_/=norm;
        sim->u[idx] += eps * (dy * Ny_ * curl[idx]);
        sim->v[idx] -= eps * (dx * Nx_ * curl[idx]);
      }

    free(curl);
}

void inject_turbulence(Simulation *sim) {
    int NX = sim->NX, M = sim->local_NY;
    double mag = sim->config.turbulence_magnitude;
    for(int j=1;j<=M;++j)
      for(int i=1;i<NX-1;++i){
        int idx = IX(i,j,NX);
        sim->u[idx] += mag * ((rand()/(double)RAND_MAX)-0.5);
        sim->v[idx] += mag * ((rand()/(double)RAND_MAX)-0.5);
      }
}

void inject_smoke(Simulation *sim,
                  int x0, int y0, int radius, double amount)
{
    int NX = sim->NX;
    double sigma = radius/2.0;
    for(int j=y0-radius; j<=y0+radius; ++j){
      for(int i=x0-radius; i<=x0+radius; ++i){
        if(i<0||i>=NX||j<1||j>sim->local_NY) continue;
        double dist = sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0));
        if(dist<=radius){
          double fac = exp(-dist*dist/(2*sigma*sigma));
          sim->density[IX(i,j,NX)] += amount * fac;
        }
      }
    }
}

void inject_heat(Simulation *sim,
                 int x0, int y0, int radius, double heat)
{
    int NX = sim->NX;
    double sigma = radius/2.0;
    for(int j=y0-radius; j<=y0+radius; ++j){
      for(int i=x0-radius; i<=x0+radius; ++i){
        if(i<0||i>=NX||j<1||j>sim->local_NY) continue;
        double dist = sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0));
        if(dist<=radius){
          double fac = exp(-dist*dist/(2*sigma*sigma));
          sim->temperature[IX(i,j,NX)] += heat * fac;
        }
      }
    }
}

// ----------------------------------------------------------------------------
// One time‐step: injection → 7 halo‐exchanges + physics → print
// ----------------------------------------------------------------------------
void simulation_step(Simulation *sim) {
    double step_start = MPI_Wtime();
    double comm_time = 0.0;
    double comp_time = 0.0;
    double temp_start;
    
    int step = (int)(sim->time / sim->dt);
    int gy   = 5;

    // 0) enforce global no-flux walls on density & temperature
    temp_start = MPI_Wtime();
    enforce_global_scalar_bc(sim, sim->density);
    enforce_global_scalar_bc(sim, sim->temperature);
    comp_time += MPI_Wtime() - temp_start;

    // 1) inject smoke/heat on rank 0
    temp_start = MPI_Wtime();
    if (sim->rank == 0) {
        if (step % sim->config.smoke_pulse_period < sim->config.smoke_pulse_duration)
            inject_smoke(sim, sim->NX/2, gy,
                         sim->config.smoke_radius,
                         sim->config.smoke_amount);
        if (step % sim->config.smoke_pulse_period == 0)
            inject_heat(sim, sim->NX/2, gy,
                        sim->config.smoke_radius, 50.0);
    }
    comp_time += MPI_Wtime() - temp_start;

    // Communication timing
    temp_start = MPI_Wtime();
    sim_exchange_ghost_rows(sim, sim->density);
    comm_time += MPI_Wtime() - temp_start;

    // 2) density advection
    temp_start = MPI_Wtime();
    double advect_start = MPI_Wtime();
    temp_start = MPI_Wtime();
    sim_exchange_ghost_rows(sim, sim->u);
    sim_exchange_ghost_rows(sim, sim->v);
    sim_exchange_ghost_rows(sim, sim->density);
    comm_time += MPI_Wtime() - temp_start;

    temp_start = MPI_Wtime();
    advect_density(sim);
    enforce_global_scalar_bc(sim, sim->density);
    sim->timing.advection_time += MPI_Wtime() - advect_start;
    comp_time += MPI_Wtime() - temp_start;

    // 3) density diffusion
    sim_exchange_ghost_rows(sim, sim->density);
    diffuse_density(sim);
    enforce_global_scalar_bc(sim, sim->density);

    // 4) temperature advection
    sim_exchange_ghost_rows(sim, sim->u);
    sim_exchange_ghost_rows(sim, sim->v);
    sim_exchange_ghost_rows(sim, sim->temperature);
    advect_temperature(sim);
    enforce_global_scalar_bc(sim, sim->temperature);

    // 5) temperature diffusion
    sim_exchange_ghost_rows(sim, sim->temperature);
    diffuse_temperature(sim);
    enforce_global_scalar_bc(sim, sim->temperature);

    // 6) velocity diffusion
    sim_exchange_ghost_rows(sim, sim->u);
    sim_exchange_ghost_rows(sim, sim->v);
    diffuse_velocity(sim);

    // 7) buoyancy, wind, turbulence, vorticity
    apply_buoyancy(sim);
    apply_wind(sim);
    if (step % sim->config.turbulence_injection_interval == 0)
        inject_turbulence(sim);
    sim_exchange_ghost_rows(sim, sim->u);
    sim_exchange_ghost_rows(sim, sim->v);
    apply_vorticity_confinement(sim);

    // 8) pressure solve
    double pressure_start = MPI_Wtime();
    temp_start = MPI_Wtime();
    sim_exchange_ghost_rows(sim, sim->u);
    sim_exchange_ghost_rows(sim, sim->v);
    sim_exchange_ghost_rows(sim, sim->pressure);
    comm_time += MPI_Wtime() - temp_start;

    temp_start = MPI_Wtime();
    solve_pressure(sim);
    sim->timing.pressure_solve_time += MPI_Wtime() - pressure_start;
    comp_time += MPI_Wtime() - temp_start;

    // 9) velocity update
    sim_exchange_ghost_rows(sim, sim->pressure);
    update_velocity(sim);

    // Update timing statistics
    sim->timing.computation_time += comp_time;
    sim->timing.communication_time += comm_time;
    sim->timing.total_time += MPI_Wtime() - step_start;
    sim->timing.step_count++;

    // advance time & print
    sim->time += sim->dt;
    if (sim->rank == 0) {
        printf("Step %d, time=%.3f\n", step, sim->time);
        printf("  Computation: %.3f ms\n", comp_time * 1000.0);
        printf("  Communication: %.3f ms\n", comm_time * 1000.0);
        printf("  Total step time: %.3f ms\n", (MPI_Wtime() - step_start) * 1000.0);
    }
}