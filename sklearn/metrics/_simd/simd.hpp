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
#endif /*SIMD_HPP*/
