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

#endif // VISUALIZATION_H
