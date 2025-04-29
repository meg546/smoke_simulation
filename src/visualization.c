#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "visualization.h"
#include "simulation.h"
#include "grid.h"

void write_vtk(Simulation *sim, const char *filename) {
    // Create the full file path in the data folder.
    char filepath[256];
    // Concatenate "data/" and the provided filename.
    snprintf(filepath, sizeof(filepath), "data/%s", filename);

    FILE *fp = fopen(filepath, "w");
    if (!fp) {
        fprintf(stderr, "Error: cannot open file %s for writing.\n", filepath);
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
