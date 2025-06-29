# MPI Parallel Smoke Simulation

A high-performance computational fluid dynamics (CFD) simulation of smoke propagation using MPI parallelization. This project implements the Navier-Stokes equations with buoyancy effects, turbulence modeling, and thermal dynamics to create realistic smoke behavior.

## Features

- **MPI Parallelization**: Efficient domain decomposition for multi-process execution
- **Physical Accuracy**: Implements Navier-Stokes equations with:
  - Velocity advection and diffusion
  - Density advection and diffusion  
  - Temperature dynamics with thermal diffusion
  - Buoyancy forces
  - Pressure projection for incompressible flow
  - Vorticity confinement for enhanced turbulence
- **Scalable Performance**: Demonstrated speedup up to 40 processes
- **Visualization**: VTK output format compatible with ParaView
- **Configurable Parameters**: Extensive simulation parameter control

## Project Structure

```
smoke_simulation_project/
├── src/                    # Source code
│   ├── main.c             # Main simulation driver
│   ├── simulation.c       # Core simulation algorithms
│   ├── simulation.h       # Simulation data structures and interfaces
│   ├── grid.c             # Grid management and utilities
│   ├── grid.h             # Grid data structures
│   ├── utils.c            # Utility functions
│   ├── utils.h            # Utility function declarations
│   ├── visualization.c    # VTK output generation
│   └── visualization.h    # Visualization interfaces
├── scripts/               # Python analysis scripts
│   └── speedup_plot.py   # Performance analysis plotting
├── results/               # Simulation outputs and plots
│   └── speedup_plot.png  # Performance analysis results
├── data/                  # Runtime data directory (VTK files)
├── examples/              # Example configurations and use cases
├── docs/                  # Additional documentation
├── Makefile              # Build configuration
└── README.md             # This file
```

## Prerequisites

### Required Software
- **MPI Implementation**: OpenMPI or MPICH
- **C Compiler**: GCC or compatible compiler with C99 support
- **Make**: Build automation
- **Python 3.7+**: For analysis scripts (optional)
- **ParaView**: For visualization (optional)

### Installing Dependencies

#### Ubuntu/Debian:
```bash
sudo apt update
sudo apt install build-essential openmpi-bin openmpi-dev python3 python3-pip
pip3 install matplotlib numpy
```

#### macOS:
```bash
brew install open-mpi python3
pip3 install matplotlib numpy
```

## Building the Project

1. **Clone the repository**:
   ```bash
   git clone https://github.com/meg546/smoke_simulation.git
   cd smoke_simulation_project
   ```

2. **Build the simulation**:
   ```bash
   make clean  # Clean previous builds
   make        # Compile the project
   ```

3. **Verify the build**:
   ```bash
   ls -la smoke_sim  # Should show the executable
   ```

## Running Simulations

### Basic Usage

**Single Process** (for testing):
```bash
mpirun -n 1 ./smoke_sim
```

**Multi-Process** (recommended):
```bash
mpirun -n 4 ./smoke_sim    # 4 processes
mpirun -n 8 ./smoke_sim    # 8 processes
mpirun -n 16 ./smoke_sim   # 16 processes
```

### Advanced Usage

**With hostfile** (for cluster execution):
```bash
mpirun -hostfile hosts.txt -n 16 ./smoke_sim
```

**With process binding**:
```bash
mpirun -n 8 --bind-to core ./smoke_sim
```

### Output

The simulation generates:
- **Console output**: Real-time performance metrics and simulation progress
- **VTK files**: `data/output_XXXXXX.vtk` files for visualization
- **Timing data**: Computation vs communication time breakdown

## Performance Analysis

### Generating Performance Plots

After running simulations with different process counts:

```bash
cd scripts
python3 speedup_plot.py
```

This generates performance analysis plots showing:
- Speedup vs number of processes
- Comparison with ideal linear speedup
- Efficiency metrics

### Benchmark Results

Current performance on a 256×256 grid:

| Processes | Speedup | Efficiency |
|-----------|---------|------------|
| 1         | 1.00    | 100%       |
| 4         | 3.54    | 88.5%      |
| 8         | 6.19    | 77.4%      |
| 16        | 9.59    | 59.9%      |
| 32        | 11.33   | 35.4%      |
| 40        | 11.84   | 29.6%      |

## Configuration Parameters

Key simulation parameters (in `simulation.c`):

### Physical Parameters
- `dt`: Time step size (default: 0.002)
- `T`: Total simulation time (default: 5.0)
- `density_diffusion`: Smoke diffusion rate (default: 0.0001)
- `velocity_viscosity`: Fluid viscosity (default: 0.0005)
- `buoyancy_beta`: Buoyancy strength (default: 0.03)

### Numerical Parameters
- `pressure_tolerance`: Pressure solver tolerance (default: 1e-5)
- `pressure_max_iter`: Maximum pressure iterations (default: 800)
- `vorticity_epsilon`: Vorticity confinement strength (default: 0.20)

### Grid Parameters (in `main.c`)
- `NX`, `NY`: Grid resolution (default: 256×256)
- `Lx`, `Ly`: Physical domain size (default: 1.0×1.0)

## Visualization

### Using ParaView

1. **Install ParaView**: Download from [paraview.org](https://www.paraview.org/)

2. **Open VTK files**:
   - Launch ParaView
   - File → Open → Navigate to `data/` directory
   - Select `output_*.vtk` files
   - Choose "Group files" if loading multiple timesteps

3. **Visualization options**:
   - **Density field**: Shows smoke concentration
   - **Velocity field**: Displays fluid motion vectors
   - **Temperature field**: Shows thermal distribution
   - **Pressure field**: Displays pressure variations

### Example Visualization Workflow
```bash
# Run simulation
mpirun -n 8 ./smoke_sim

# Files will be in data/ directory
ls data/output_*.vtk

# Open in ParaView for visualization
paraview data/output_000000.vtk
```

## Algorithm Details

### Navier-Stokes Implementation

The simulation uses a projection method with the following steps:

1. **Advection**: Semi-Lagrangian method with wall reflection
2. **Diffusion**: Explicit finite difference with boundary conditions
3. **External Forces**: Buoyancy, wind, and turbulence injection
4. **Pressure Projection**: Gauss-Seidel iteration with Neumann boundaries
5. **Velocity Update**: Pressure gradient correction

### MPI Parallelization Strategy

- **Domain Decomposition**: Horizontal strips across processes
- **Ghost Rows**: Two-row halo exchange for stencil operations
- **Load Balancing**: Automatic distribution of remainder rows
- **Communication Pattern**: Nearest-neighbor exchange with `MPI_Sendrecv`

### Debug Mode

Compile with debug symbols:
```bash
make clean
CFLAGS="-Wall -g -O0 -std=c99" make
```

Run with debugger:
```bash
mpirun -n 4 gdb ./smoke_sim
```

