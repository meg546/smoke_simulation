#!/bin/bash
# Example script for running performance benchmarks

echo "Smoke Simulation Performance Benchmark"
echo "======================================"

# Check if executable exists
if [ ! -f "./smoke_sim" ]; then
    echo "Error: smoke_sim executable not found. Please run 'make' first."
    exit 1
fi

# Create data directory if it doesn't exist
mkdir -p data results

# Array of process counts to test
PROCESS_COUNTS=(1 2 4 8 16)

echo "Running benchmarks with different process counts..."

for n in "${PROCESS_COUNTS[@]}"; do
    echo ""
    echo "Running with $n process(es)..."
    echo "------------------------------"
    
    # Clear previous data
    rm -f data/output_*.vtk
    
    # Run simulation and capture timing
    start_time=$(date +%s.%N)
    mpirun -n $n ./smoke_sim > "results/timing_${n}proc.log" 2>&1
    end_time=$(date +%s.%N)
    
    # Calculate elapsed time
    elapsed=$(echo "$end_time - $start_time" | bc -l)
    echo "Completed in ${elapsed} seconds"
    
    # Count output files
    vtk_count=$(ls data/output_*.vtk 2>/dev/null | wc -l)
    echo "Generated $vtk_count VTK files"
done

echo ""
echo "Benchmark complete!"
echo "Results saved in results/ directory"
echo "Run 'python3 scripts/speedup_plot.py' to generate performance plots"
