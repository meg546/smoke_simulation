#ifndef GRID_H
#define GRID_H

typedef struct {
    int NX;
    int NY;
    double dx;
    double dy;
    double *x_coords;
    double *y_coords;
} Grid;

/**
 * Creates a grid covering a domain of size Lx by Ly with NX and NY grid points.
 * Returns a pointer to the allocated Grid structure, or NULL if allocation fails.
 */
Grid* create_grid(int NX, int NY, double Lx, double Ly);

void free_grid(Grid *grid);

void print_grid(const Grid *grid);

#endif