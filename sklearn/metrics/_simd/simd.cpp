#include "simd.hpp"

template float simd_manhattan_dist(
    const float* HWY_RESTRICT x,
    const float* HWY_RESTRICT y,
    const size_t size
);
template double simd_manhattan_dist(
    const double* HWY_RESTRICT x,
    const double* HWY_RESTRICT y,
    const size_t size
);
