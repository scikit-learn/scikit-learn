#ifndef SIMD_HPP
#define SIMD_HPP
#include <hwy/base.h>

namespace manhattan{

    template <typename Type>
    HWY_DLLEXPORT Type simd_manhattan_dist(
        const Type* x,
        const Type* y,
        const size_t size
    );
}

namespace manhattan{
    namespace N_AVX2{
        float manhattan_dist_float(
            const float* x,
            const float* y,
            const size_t size
        );
        double manhattan_dist_double(
            const double* x,
            const double* y,
            const size_t size
        );
    }
}
#endif /*SIMD_HPP*/
