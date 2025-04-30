# Compiler and flags
CC = mpicc
CFLAGS = -Wall -O3 -std=c99

# Source files (update if you add more modules)
SRCS = src/main.c src/grid.c src/simulation.c src/utils.c src/visualization.c src/mpi_utils.c

# Object files
OBJS = $(SRCS:.c=.o)

# Target executable
TARGET = smoke_simulation

# Default rule
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule: Remove object files, the executable, and all files in the data folder.
ifeq ($(OS),Windows_NT)
    RM = del /Q
    # For Windows, list files using backslashes and wildcard pattern.
    CLEAN_FILES = src\main.o src\grid.o src\simulation.o src\utils.o src\visualization.o src\mpi_utils.o $(TARGET).exe output.vtk data\*
else
    RM = rm -f
    CLEAN_FILES = $(OBJS) $(TARGET) output.vtk data/*
endif

clean:
	$(RM) $(CLEAN_FILES)

# Run rule for convenience
run: $(TARGET)
	mpirun -np 4 ./$(TARGET)

.PHONY: all clean run
