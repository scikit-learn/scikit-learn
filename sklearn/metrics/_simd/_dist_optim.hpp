/*
Only provide a body upon inclusion if we are being included in _dist_metrics.pyx
otherwise we remain empty so as to bypass Cython's forced inclusion of
this file due to cimporting _dist_metrics
*/
#if defined(DIST_METRICS)
#include "simd.hpp"

/* If building with SIMD support, include the compiled library code */
#if WITH_SIMD != 1
#include <cstddef>
#include "hwy/base.h"

/* Else, we provide trivial functions for compilation */
namespace manhattan{
    HWY_DLLEXPORT float simd_manhattan_dist_f32(
        const float* x,
        const float* y,
        const size_t size
    ){return -1;}
    HWY_DLLEXPORT double simd_manhattan_dist_f64(
        const double* x,
        const double* y,
        const size_t size
    ){return -1;}
    HWY_DLLEXPORT func32_t get_simd_manhattan_f32(){
        return nullptr;
    }
    HWY_DLLEXPORT func64_t get_simd_manhattan_f64(){
        return nullptr;
    }
}
#endif
#else
/* Empty body */
#endif
