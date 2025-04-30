#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include "simulation.h"

/**
 * Write the simulation fields to a legacy VTK file.
 * The file contains a structured grid with:
 *  - POINTS: coordinates of each grid point (z set to 0)
 *  - POINT_DATA: includes scalar density and vector velocity.
 */
void write_vtk(Simulation *sim, const char *filename);
void write_vtk_global(int NX, int global_NY,
    double dx, double dy,
    const double *global_density,
    const char *filename);

#endif // VISUALIZATION_H
