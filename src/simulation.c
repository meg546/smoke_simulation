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
    c.dt                             = 0.005;
    c.T                              = 5.0;
    c.density_diffusion              = 0.0001;
    c.velocity_viscosity             = 0.005;
    c.pressure_tolerance             = 1e-6;
    c.pressure_max_iter              = 1000;
    c.buoyancy_beta                  = 0.08;
    c.ambient_density                = 0.01;
    c.wind_u                         = 4.5;
    c.wind_v                         = 3.5;
    c.turbulence_magnitude           = 0.02;
    c.turbulence_injection_interval  = 2;
    c.vorticity_epsilon              = 0.2;
    c.smoke_radius                   = 8;
    c.smoke_amount                   = 0.4;
    c.smoke_pulse_period             = 10;
    c.smoke_pulse_duration           = 5;
    c.velocity_clamp_min             = -5.0;
    c.velocity_clamp_max             = 5.0;
    c.thermal_diffusivity            = 0.0001;
    c.ambient_temperature            = 300.0;
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

    return sim;
}

void free_simulation(Simulation *sim) {
    if (!sim) return;
    free(sim->density);
    free(sim->u);
    free(sim->v);
    free(sim->pressure);
    free(sim->temperature);
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

    // preserve ghost‐rows so we only overwrite interior+ghosts in memcpy:
    memcpy(newD, old, NX*(M+2)*sizeof(double));

    for (int j = 1; j <= M; ++j) {
        double y = (sim->j0 + j - 1) * dy;
        for (int i = 0; i < NX; ++i) {
            int idx = IX(i,j,NX);
            double x = i*dx;
            double u = old[idx], v = sim->v[idx];
            double xs = fmax(0.0, fmin(x - u*dt, (NX-1)*dx));
            double ys = fmax(0.0, fmin(y - v*dt, (sim->NY-1)*dy));

            int i0  = (int)floor(xs/dx),
                j0g = (int)floor(ys/dy);
            int i1  = (i0+1 < NX ? i0+1 : NX-1),
                j1g = (j0g+1 < sim->NY ? j0g+1 : sim->NY-1);

            // map global j into local index [0..M+1]
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
    for(int it=0; it<maxit; ++it){
      double maxchg=0;
      for(int j=1;j<=M;++j){
        for(int i=1;i<NX-1;++i){
          int idx = IX(i,j,NX);
          double div = (sim->u[IX(i+1,j,NX)] - sim->u[IX(i-1,j,NX)])/(2*dx)
                     + (sim->v[IX(i,j+1,NX)] - sim->v[IX(i,j-1,NX)])/(2*dy);
          double val = ( sim->pressure[IX(i-1,j,NX)]
                       + sim->pressure[IX(i+1,j,NX)]
                       + sim->pressure[IX(i,j-1,NX)]
                       + sim->pressure[IX(i,j+1,NX)]
                       - (dx*dx/dt)*div ) * 0.25;
          maxchg = fmax(maxchg, fabs(val - sim->pressure[idx]));
          p2[idx] = val;
        }
      }
      memcpy(&sim->pressure[IX(1,1,NX)],
             &p2[IX(1,1,NX)],
             (NX-2)*M*sizeof(double));
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

    for(int j=1;j<=M;++j){
      for(int i=1;i<NX-1;++i){
        int idx = IX(i,j,NX);
        double dpdx = (sim->pressure[IX(i+1,j,NX)] - sim->pressure[IX(i-1,j,NX)])/(2*dx);
        double dpdy = (sim->pressure[IX(i,j+1,NX)] - sim->pressure[IX(i,j-1,NX)])/(2*dy);
        u2[idx] = clamp(sim->u[idx] - dt*dpdx, vmin, vmax);
        v2[idx] = clamp(sim->v[idx] - dt*dpdy, vmin, vmax);
      }
    }
    // no‐slip walls
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
    int step = (int)(sim->time / sim->dt);
    int gy   = 5;

    // 0) enforce global no-flux walls on density & temperature
    enforce_global_scalar_bc(sim, sim->density);
    enforce_global_scalar_bc(sim, sim->temperature);

    // 1) inject smoke/heat on rank 0
    if (sim->rank == 0) {
        if (step % sim->config.smoke_pulse_period < sim->config.smoke_pulse_duration)
            inject_smoke(sim, sim->NX/2, gy,
                         sim->config.smoke_radius,
                         sim->config.smoke_amount);
        if (step % sim->config.smoke_pulse_period == 0)
            inject_heat(sim, sim->NX/2, gy,
                        sim->config.smoke_radius, 50.0);
        // push into rank-1’s ghost immediately
        
    }
    sim_exchange_ghost_rows(sim, sim->density);

    // 2) density advection
    sim_exchange_ghost_rows(sim, sim->u);
    sim_exchange_ghost_rows(sim, sim->v);
    sim_exchange_ghost_rows(sim, sim->density);
    advect_density(sim);
    enforce_global_scalar_bc(sim, sim->density);

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
    sim_exchange_ghost_rows(sim, sim->u);
    sim_exchange_ghost_rows(sim, sim->v);
    sim_exchange_ghost_rows(sim, sim->pressure);
    solve_pressure(sim);

    // 9) velocity update
    sim_exchange_ghost_rows(sim, sim->pressure);
    update_velocity(sim);

    // advance time & print
    sim->time += sim->dt;
    if (sim->rank == 0)
        printf("Completed step %d, time=%.3f\n", step, sim->time);
}