#include <stdio.h>
#include <stdlib.h>
#include "grid.h"

Grid* create_grid(int NX, int NY, double Lx, double Ly)
{
    Grid* grid = (Grid*)malloc(sizeof(Grid));
    if (!grid)
    {
        fprintf(stderr, "Error allocating memory for Grid");
        return NULL;
    }

    grid->NX = NX;
    grid->NY = NY;
    grid->dx = Lx / (NX - 1);
    grid->dy = Ly / (NY - 1);

    grid->x_coords = (double*)malloc(NX * sizeof(double));
    grid->y_coords = (double*)malloc(NY * sizeof(double));
    if (!grid->x_coords || !grid->y_coords) {
        fprintf(stderr, "Error allocating memory for grid coordinates\n");
        free(grid->x_coords);
        free(grid->y_coords);
        free(grid);
        return NULL;
    } 

    for (int i = 0; i < NX; i++)
    {
        grid->x_coords[i] = i * grid->dx;
    }
    for (int j = 0; j < NY; j++) {
        grid->y_coords[j] = j * grid->dy;
    }


    return grid;
}

void print_grid(const Grid *grid) {
    if (!grid) {
        printf("Grid is NULL.\n");
        return;
    }
    printf("Grid Dimensions: %d x %d\n", grid->NX, grid->NY);
    printf("Grid Spacing: dx = %.4f, dy = %.4f\n", grid->dx, grid->dy);
    printf("Grid Coordinates:\n");
    for (int j = 0; j < grid->NY; j++) {
        for (int i = 0; i < grid->NX; i++) {
            printf("(%.2f, %.2f) ", grid->x_coords[i], grid->y_coords[j]);
        }
        printf("\n");
    }
}

void free_grid(Grid *grid) {
    if (grid) {
        free(grid->x_coords);
        free(grid->y_coords);
        free(grid);
    }
}

