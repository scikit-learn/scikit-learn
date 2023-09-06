/*
Only provide a body upon inclusion if we are being included in _dist_metrics.pyx
otherwise we remain empty so as to bypass Cython's forced inclusion of
this file due to cimporting _dist_metrics
*/
#if defined(DIST_METRICS)

/* If built with SIMD support, include the compiled library code */
#if WITH_SIMD == 1
#include "simd.hpp"
#else
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
}

/*
In case of a runtime machine without AVX, we need to
provide alternative scalar implementations.
*/
template <typename Type>
Type simd_manhattan_dist_scalar(
    const Type* x,
    const Type* y,
    const size_t size
){
    double scalar_sum = 0;

    for(std::size_t idx = 0; idx < size; ++idx) {
        scalar_sum += fabs(x[idx] - y[idx]);
    }
    return (Type) scalar_sum;
}
#endif
#else
/* Empty body */
#endif
