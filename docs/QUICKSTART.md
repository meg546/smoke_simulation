# Quick Start Guide

## 5-Minute Setup

### 1. Install Dependencies
```bash
# Ubuntu/Debian
sudo apt install build-essential openmpi-bin openmpi-dev

# macOS
brew install open-mpi
```

### 2. Build and Run
```bash
# Clone and build
git clone https://github.com/meg546/smoke_simulation.git
cd smoke_simulation_project
make

# Quick test
make test

# Full simulation
make run-4
```

### 3. View Results
```bash
# List output files
ls data/output_*.vtk

# Open in ParaView (if installed)
paraview data/output_000000.vtk
```

## Common Commands

```bash
# Build project
make

# Run with different process counts
make run-4    # 4 processes
make run-8    # 8 processes

# Performance analysis
make analyze

# Clean up
make clean

# Help
make help
```

## Troubleshooting

**"mpirun not found"**
```bash
# Check MPI installation
which mpirun
mpicc --version
```

**Build errors**
```bash
# Clean and rebuild
make clean
make
```

**No output files**
```bash
# Check if data directory exists
ls -la data/
# Should contain output_*.vtk files after running
```

## Next Steps

1. **Modify parameters**: See `examples/config_examples.c`
2. **Performance testing**: Run `examples/run_benchmark.sh`
3. **Visualization**: Install ParaView for advanced visualization
4. **Read full documentation**: See main `README.md`
