#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "visualization.h"
#include "simulation.h"
#include "grid.h"

void write_vtk(Simulation *sim, const char *filename) {
    // Create full filepath with data/ prefix
    char filepath[256];
    snprintf(filepath, sizeof(filepath), "data/%s", filename);

    // Try to create data directory if it doesn't exist
    if (system("mkdir -p data") != 0) {
        fprintf(stderr, "Warning: Failed to create data directory\n");
    }

    FILE *fp = fopen(filepath, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s for writing.\n", filepath);
        return;
    }
    
    int NX = sim->NX;
    int NY = sim->NY;
    int num_points = NX * NY;
    double dx = sim->grid->dx;
    double dy = sim->grid->dy;
    
    // Write VTK header information.
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Smoke Simulation Data\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d 1\n", NX, NY);
    fprintf(fp, "POINTS %d float\n", num_points);
    
    // Write point coordinates.
    // Assumes grid is organized in row-major order: for each j (y), for each i (x).
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            double x = i * dx;
            double y = j * dy;
            fprintf(fp, "%f %f %f\n", x, y, 0.0);
        }
    }
    
    // Write point data header.
    fprintf(fp, "\nPOINT_DATA %d\n", num_points);
    
    // Write scalar field: density.
    fprintf(fp, "SCALARS density float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            int idx = IX(i, j, NX);
            fprintf(fp, "%f\n", sim->density[idx]);
        }
    }
    
    // Write vector field: velocity (with z component = 0).
    fprintf(fp, "\nVECTORS velocity float\n");
    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            int idx = IX(i, j, NX);
            fprintf(fp, "%f %f %f\n", sim->u[idx], sim->v[idx], 0.0);
        }
    }
    
    fclose(fp);
}


void write_vtk_global(int NX, int global_NY,
    double dx, double dy,
    const double *global_density,
    const char *filename)
{
    // Create full filepath with data/ prefix
    char filepath[256];
    snprintf(filepath, sizeof(filepath), "data/%s", filename);

    FILE *f = fopen(filepath, "w");
    if (!f) 
    {
        // Create data directory if it doesn't exist
        if (system("mkdir -p data") != 0) {
            fprintf(stderr, "Warning: Failed to create data directory\n");
        }
        // Try opening file again
        f = fopen(filepath, "w");
        if (!f) {
            fprintf(stderr, "Error: Cannot open file %s for writing.\n", filepath);
            return;
        }
    }

    fprintf(f,
    "# vtk DataFile Version 3.0\n"
    "Smoke density (gathered)\n"
    "ASCII\n"
    "DATASET STRUCTURED_POINTS\n"
    "DIMENSIONS %d %d 1\n"
    "ORIGIN 0 0 0\n"
    "SPACING %g %g 1.0\n"
    "POINT_DATA %d\n"
    "SCALARS density double 1\n"
    "LOOKUP_TABLE default\n",
    NX, global_NY,
    dx, dy,
    NX*global_NY
    );

    for (int j = 0; j < global_NY; j++) 
    {
        for (int i = 0; i < NX; i++) 
        {
            fprintf(f, "%g\n", global_density[j*NX + i]);
        }
    }
    
    fclose(f);
}
