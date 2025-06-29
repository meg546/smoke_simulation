# Makefile for MPI-parallel smoke simulation

# Compiler and flags
CC       = mpicc
CFLAGS   = -Wall -O3 -std=c99
LDFLAGS  = -lm

# Source directory
SRC_DIR  = src

# Source files
SRCS     = $(SRC_DIR)/main.c \
           $(SRC_DIR)/grid.c \
           $(SRC_DIR)/simulation.c \
           $(SRC_DIR)/utils.c \
           $(SRC_DIR)/visualization.c

# Object files (one .o per .c)
OBJS     = $(SRCS:.c=.o)

# Target executable
TARGET   = smoke_sim

# Default rule
all: $(TARGET)

# Link step
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile step
$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule: removes object files, executable, vtk outputs, and data/*
clean:
	@echo "Cleaning..."
	-rm -f $(OBJS) $(TARGET)
	-rm -f out_r*_step*.vtk
	-rm -rf data/*

# Clean everything including results
clean-all: clean
	@echo "Cleaning all outputs..."
	-rm -rf results/*.png results/*.dat

# Create necessary directories
setup:
	@echo "Setting up project directories..."
	@mkdir -p data results docs examples

# Run performance analysis
analyze:
	@echo "Running performance analysis..."
	cd scripts && python3 speedup_plot.py

# Quick test with single process
test:
	@echo "Running single-process test..."
	mpirun -n 1 ./$(TARGET)

# Example multi-process runs
run-4:
	@echo "Running with 4 processes..."
	mpirun -n 4 ./$(TARGET)

run-8:
	@echo "Running with 8 processes..."
	mpirun -n 8 ./$(TARGET)

# Help target
help:
	@echo "Available targets:"
	@echo "  all        - Build the simulation (default)"
	@echo "  clean      - Remove build files and data"
	@echo "  clean-all  - Remove all build files and results"
	@echo "  setup      - Create necessary directories"
	@echo "  test       - Run single-process test"
	@echo "  run-4      - Run with 4 processes"
	@echo "  run-8      - Run with 8 processes"
	@echo "  analyze    - Generate performance plots"
	@echo "  help       - Show this help message"

.PHONY: all clean clean-all setup analyze test run-4 run-8 help
