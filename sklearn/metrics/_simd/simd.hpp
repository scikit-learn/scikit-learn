#ifndef SIMD_HPP
#define SIMD_HPP
#include <hwy/base.h>

namespace manhattan{
    using func32_t = float (*)(const float *x, const float *y, size_t size);
    using func64_t = double (*)(const double *x, const double *y, size_t size);

    HWY_DLLEXPORT float simd_manhattan_dist_f32(
        const float* x,
        const float* y,
        const size_t size
    );
    HWY_DLLEXPORT double simd_manhattan_dist_f64(
        const double* x,
        const double* y,
        const size_t size
    );
    HWY_DLLEXPORT func32_t get_simd_manhattan_f32();
    HWY_DLLEXPORT func64_t get_simd_manhattan_f64();
}
#endif /*SIMD_HPP*/
