# Smoke Simulation

A parallel computational fluid dynamics simulation for smoke dispersion using MPI for distributed computing.

## Overview

This project simulates smoke dynamics in a 2D fluid environment using numerical methods to solve the Navier-Stokes equations. The simulation includes:

- **Fluid Dynamics**: Velocity field computation with viscosity effects
- **Density Transport**: Smoke/gas concentration advection and diffusion
- **Buoyancy Forces**: Temperature-driven vertical motion
- **Pressure Projection**: Incompressible flow enforcement
- **Vorticity Confinement**: Enhanced turbulence and swirl preservation
- **Parallel Computing**: MPI-based domain decomposition for high-performance computing

## Features

### Physical Phenomena
- **Advection**: Transport of smoke by fluid velocity
- **Diffusion**: Spreading of smoke concentration over time
- **Buoyancy**: Temperature-driven upward motion of hot gases
- **Turbulence**: Chaotic flow patterns and vortex dynamics
- **Wind Effects**: External wind forces on the simulation
- **Vorticity Confinement**: Artificial enhancement of rotational flow structures

### Computational Features
- **MPI Parallelization**: Distributed computing across multiple processors
- **Domain Decomposition**: Grid partitioning for parallel processing
- **Ghost Cell Communication**: Boundary data exchange between processes
- **VTK Output**: Industry-standard visualization format
- **Configurable Parameters**: Adjustable physical and numerical parameters

## Requirements

### Dependencies
- **MPI Implementation**: OpenMPI, MPICH, or Intel MPI
- **C Compiler**: GCC, Clang, or Intel C Compiler with C99 support
- **Math Library**: Standard C math library (-lm)

### System Requirements
- Unix-like system (Linux, macOS) or Windows with MPI support
- Multiple CPU cores for parallel execution (recommended)
- Sufficient memory for grid storage (depends on grid size)

## Installation

### 1. Clone the Repository
```bash
git clone <repository-url>
cd smoke_simulation
```

### 2. Install MPI
#### On Ubuntu/Debian:
```bash
sudo apt-get install libopenmpi-dev openmpi-bin
```

#### On macOS (with Homebrew):
```bash
brew install open-mpi
```

#### On Windows:
Install Microsoft MPI or Intel MPI runtime and SDK.

### 3. Build the Project
```bash
make
```

This will compile the source files and create the `smoke_simulation` executable.

## Usage

### Basic Execution
```bash
# Run with 4 MPI processes
mpirun -np 4 ./smoke_simulation
```

### Alternative MPI Launchers
```bash
# Intel MPI
mpiexec -n 4 ./smoke_simulation

# SLURM (on HPC clusters)
srun -n 4 ./smoke_simulation
```

### Build Targets
```bash
make all      # Build the executable (default)
make clean    # Remove object files and executables
make run      # Build and run with 4 processes
```

## Configuration

The simulation parameters can be modified in `src/main.c`:

### Grid Parameters
```c
int NX = 100;           // Grid cells in x-direction
int NY = 100;           // Grid cells in y-direction
double Lx = 1.0;        // Domain length in X
double Ly = 1.0;        // Domain length in Y
```

### Time Parameters
```c
double dt = 0.005;      // Time step size
double T = 5.0;         // Total simulation time
int output_interval = 1; // Output frequency
```

### Physical Parameters (in `simulation.h`)
- `density_diffusion`: Controls smoke spreading rate
- `velocity_viscosity`: Fluid viscosity (resistance to flow)
- `buoyancy_beta`: Strength of buoyancy forces
- `wind_u`, `wind_v`: External wind velocity
- `turbulence_magnitude`: Turbulence injection strength
- `vorticity_epsilon`: Vorticity confinement strength

## Output

The simulation generates VTK files for visualization:
- **File Format**: Legacy VTK format (.vtk)
- **Naming Convention**: `output_XXXX.vtk` (where XXXX is the time step)
- **Data Fields**:
  - Density (scalar field)
  - Velocity (vector field)
  - Grid coordinates

### Visualization Tools
- **ParaView**: Professional scientific visualization (recommended)
- **VisIt**: Open-source visualization platform
- **VTK-based viewers**: Any software supporting VTK format

## Project Structure

```
smoke_simulation/
├── src/                    # Source files
│   ├── main.c             # Main program and simulation loop
│   ├── grid.c             # Grid creation and management
│   ├── simulation.c       # Core simulation physics
│   ├── mpi_utils.c        # MPI communication routines
│   ├── utils.c            # Utility functions
│   └── visualization.c    # VTK output generation
├── include/               # Header files
│   ├── grid.h
│   ├── simulation.h
│   ├── mpi_utils.h
│   ├── utils.h
│   └── visualization.h
├── build/                 # Build artifacts (object files, executables)
├── data/                  # Output data files
├── Makefile              # Build configuration
├── .gitignore            # Git ignore rules
└── README.md             # This file
```

## Algorithm Details

### Numerical Methods
1. **Semi-Lagrangian Advection**: Stable transport of quantities
2. **Explicit Diffusion**: Forward Euler method for diffusion
3. **Pressure Projection**: Iterative solver for incompressible flow
4. **Finite Difference**: Spatial derivatives on structured grid

### MPI Parallelization
1. **Domain Decomposition**: 2D grid partitioned among processes
2. **Ghost Cell Exchange**: Communication of boundary data
3. **Collective Operations**: Global data gathering for output
4. **Load Balancing**: Even distribution of computational work

### Simulation Steps (per time step)
1. Exchange ghost cells between MPI processes
2. Apply external forces (gravity, wind, turbulence)
3. Advect velocity field
4. Apply viscosity (velocity diffusion)
5. Solve pressure projection (incompressible flow)
6. Advect density field
7. Apply density diffusion
8. Output results (if scheduled)

## Performance Notes

### Scaling Recommendations
- **Grid Size**: Larger grids benefit more from parallelization
- **Process Count**: Optimal number depends on grid size and system
- **Memory Usage**: Approximately `8 * NX * NY * bytes_per_double` per process

### Optimization Tips
- Use compiler optimization flags (`-O3`)
- Ensure good load balancing (square grid decomposition)
- Monitor MPI communication overhead
- Consider cache-friendly memory access patterns

## Troubleshooting

### Common Issues
1. **MPI Not Found**: Ensure MPI is properly installed and in PATH
2. **Compilation Errors**: Check compiler version and flags
3. **Runtime Crashes**: Verify grid size fits in available memory
4. **Poor Performance**: Check process count vs. problem size ratio

### Debug Mode
Compile with debug flags for troubleshooting:
```bash
make CFLAGS="-Wall -g -DDEBUG"
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request


## Contact

meg546@nau.edu
https://www.linkedin.com/in/matthewgardne/

