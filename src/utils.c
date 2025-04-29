#include "utils.h"

// Helper functions to compute min and max over an array.
double min_array(const double *arr, int size) {
    double min_val = arr[0];
    for (int i = 1; i < size; i++) {
        if (arr[i] < min_val)
            min_val = arr[i];
    }
    return min_val;
}

double max_array(const double *arr, int size) {
    double max_val = arr[0];
    for (int i = 1; i < size; i++) {
        if (arr[i] > max_val)
            max_val = arr[i];
    }
    return max_val;
}

// Clamping helper function.
double clamp(double val, double lower, double upper) {
    if (val < lower)
        return lower;
    else if (val > upper)
        return upper;
    else
        return val;
}