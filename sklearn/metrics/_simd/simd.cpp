#include "simd.hpp"

template float simd_manhattan_dist(
    const float* x,
    const float* y,
    const size_t size
);
template double simd_manhattan_dist(
    const double* x,
    const double* y,
    const size_t size
);
