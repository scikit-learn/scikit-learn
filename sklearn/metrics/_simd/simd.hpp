#ifndef SIMD_HPP
#define SIMD_HPP
#include <hwy/base.h>

namespace manhattan{

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
}
#endif /*SIMD_HPP*/
