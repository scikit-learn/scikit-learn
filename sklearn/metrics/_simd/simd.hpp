#ifndef SIMD_HPP
#define SIMD_HPP
#include <hwy/base.h>

namespace manhattan{

    HWY_DLLEXPORT float _simd_manhattan_dist(
        const float* x,
        const float* y,
        const size_t size
    );
}
#endif /*SIMD_HPP*/
