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

.PHONY: all clean
