# Technical Documentation

## Algorithm Overview

The smoke simulation implements a projection method for solving the incompressible Navier-Stokes equations with additional physics for realistic smoke behavior.

### Governing Equations

**Navier-Stokes (Momentum):**
```
∂u/∂t + (u·∇)u = -∇p + ν∇²u + f
```

**Continuity (Incompressibility):**
```
∇·u = 0
```

**Density Transport:**
```
∂ρ/∂t + (u·∇)ρ = D∇²ρ + S_ρ
```

**Temperature Transport:**
```
∂T/∂t + (u·∇)T = α∇²T + S_T
```

Where:
- `u`: velocity field
- `p`: pressure
- `ρ`: density (smoke concentration)
- `T`: temperature
- `ν`: kinematic viscosity
- `D`: density diffusion coefficient
- `α`: thermal diffusivity
- `f`: external forces (buoyancy, wind)
- `S_ρ`, `S_T`: source terms

### Numerical Method

**Time Integration:** Explicit Euler with operator splitting
**Spatial Discretization:** Finite differences on regular grid
**Advection:** Semi-Lagrangian with reflection boundaries
**Pressure Solve:** Gauss-Seidel iteration

### Algorithm Steps (per timestep)

1. **Advection** (semi-Lagrangian):
   - Backtrack particle positions
   - Interpolate field values
   - Apply reflection at boundaries

2. **Diffusion** (explicit finite differences):
   - Apply 5-point Laplacian stencil
   - Forward Euler integration
   - Enforce boundary conditions

3. **External Forces**:
   - Buoyancy: `f_y = β(T - T_amb)`
   - Wind: constant velocity addition
   - Turbulence: random perturbations

4. **Vorticity Confinement**:
   - Compute vorticity: `ω = ∇×u`
   - Calculate confinement force
   - Add to velocity field

5. **Pressure Projection**:
   - Solve Poisson equation: `∇²p = ∇·u/Δt`
   - Update velocity: `u_new = u - Δt∇p`

## MPI Parallelization

### Domain Decomposition
- **Strategy**: Horizontal strips (1D decomposition in y-direction)
- **Load balancing**: Distribute remainder rows evenly
- **Ghost cells**: 2-row halo for finite difference stencils

### Communication Pattern
```
Process layout (4 processes):
┌─────────────┐ ← rank 3 (top)
├─────────────┤
├─────────────┤ ← rank 2
├─────────────┤
├─────────────┤ ← rank 1
├─────────────┤
└─────────────┘ ← rank 0 (bottom)
```

### Ghost Row Exchange
```c
// Send bottom boundary to rank+1, receive from rank+1
// Send top boundary to rank-1, receive from rank-1
MPI_Sendrecv(send_buffer, count, MPI_DOUBLE, neighbor, tag,
             recv_buffer, count, MPI_DOUBLE, neighbor, tag,
             MPI_COMM_WORLD, &status);
```

### Boundary Conditions
- **Velocity**: No-slip walls (u=v=0)
- **Scalar fields**: No-flux (Neumann) boundaries
- **Global walls**: Only handled by boundary processes
- **MPI interfaces**: Handled by ghost row exchange

## Data Structures

### Grid Structure
```c
typedef struct {
    int NX, NY;          // Grid dimensions
    double Lx, Ly;       // Physical domain size
    double dx, dy;       // Grid spacing
} Grid;
```

### Simulation State
```c
typedef struct {
    Grid *grid;
    double *density;     // Smoke concentration
    double *u, *v;       // Velocity components
    double *pressure;    // Pressure field
    double *temperature; // Temperature field
    
    int rank, nprocs;    // MPI info
    int local_NY, j0;    // Local domain info
    // ... configuration and timing data
} Simulation;
```

## Performance Characteristics

### Computational Complexity
- **Per timestep**: O(N²) for N×N grid
- **Memory usage**: O(N²/P) per process for P processes
- **Communication**: O(N) boundary exchange per process

### Scalability Bottlenecks
1. **Communication overhead**: Increases with process count
2. **Load imbalance**: Minimal for regular grids
3. **Memory bandwidth**: Limits single-node performance
4. **Algorithm**: Explicit time stepping limits parallelism

### Optimization Strategies
- **Cache optimization**: Iterate in memory order
- **Vectorization**: Use compiler auto-vectorization
- **Communication overlap**: Overlap computation with MPI
- **Load balancing**: Even distribution of grid points

## Numerical Stability

### CFL Condition (Advection)
```
Δt ≤ min(Δx/|u_max|, Δy/|v_max|)
```

### Diffusion Stability
```
Δt ≤ min(Δx²/(2D), Δy²/(2D))  # For diffusion coefficient D
```

### Recommended Parameters
- Grid: 256×256 for good detail/performance balance
- Time step: 0.002-0.003 seconds
- Diffusion coefficients: O(10⁻⁴) for stability

## File Formats

### VTK Output Structure
```
# vtk DataFile Version 3.0
Smoke Simulation Output
ASCII
DATASET STRUCTURED_POINTS
DIMENSIONS nx ny 1
SPACING dx dy 1.0
ORIGIN 0.0 0.0 0.0

POINT_DATA npts
SCALARS density double 1
LOOKUP_TABLE default
[density values...]

VECTORS velocity double
[velocity values...]
```

### Data Organization
- Each timestep → separate VTK file
- Fields: density, velocity, pressure, temperature
- Compatible with ParaView, VisIt, other VTK viewers
