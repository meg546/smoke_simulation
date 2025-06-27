#ifndef UTILS_H
#define UTILS_H

// Helper functions to compute min and max over an array.
double min_array(const double *arr, int size);
double max_array(const double *arr, int size);
// Clamping helper function.
double clamp(double val, double lower, double upper);

#endif
