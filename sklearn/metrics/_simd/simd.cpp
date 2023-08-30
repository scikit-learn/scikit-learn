#include "simd.hpp"
#include <cmath>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "simd.cpp"
#include "hwy/foreach_target.h"
#include "hwy/highway.h"
HWY_BEFORE_NAMESPACE();  // at file scope

namespace manhattan {

    namespace HWY_NAMESPACE {
        namespace hn = hwy::HWY_NAMESPACE;


        template <typename Type>
        inline Type manhattan_dist(
            const Type* x,
            const Type* y,
            const size_t size
        ) {
            printf("DEBUG *** %s\n", hwy::TargetName(HWY_TARGET));
            const hn::ScalableTag<Type> d;
            auto simd_sum_1 = hn::Zero(d);
            auto simd_sum_2 = hn::Zero(d);

            auto lane_size = hn::Lanes(d);
            size_t loop_iter = lane_size * 2;
            size_t vec_size = size - size % loop_iter;
            size_t vec_remainder_size = size - size % lane_size;

            for (size_t i = 0; i < vec_size; i += loop_iter) {
                const auto simd_x_1 = hn::LoadU(d, x + i);
                const auto simd_y_1 = hn::LoadU(d, y + i);
                simd_sum_1 += hn::AbsDiff(simd_x_1, simd_y_1);

                const auto simd_x_2 = hn::LoadU(d, x + i + lane_size);
                const auto simd_y_2 = hn::LoadU(d, y + i + lane_size);
                simd_sum_2 += hn::AbsDiff(simd_x_2, simd_y_2);
            }
            for (size_t i = vec_size; i < vec_remainder_size; i += loop_iter) {
                const auto simd_x_1 = hn::LoadU(d, x + i);
                const auto simd_y_1 = hn::LoadU(d, y + i);
                simd_sum_1 += hn::AbsDiff(simd_x_1, simd_y_1);
            }
            simd_sum_1 += simd_sum_2;
            Type scalar_sum = hn::ReduceSum(d, simd_sum_1);
            for (size_t i = vec_remainder_size; i < size; i += 1) {
                scalar_sum += fabs(x[i] - y[i]);
            }
            return scalar_sum;
        }
        float manhattan_dist_float(
            const float* x,
            const float* y,
            const size_t size
        ) {
            return manhattan_dist<float>(x, y, size);
        }
        double manhattan_dist_double(
            const double* x,
            const double* y,
            const size_t size
        ) {
            return manhattan_dist<double>(x, y, size);
        }
    }
}
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace manhattan {

    HWY_EXPORT(manhattan_dist_float);
    HWY_EXPORT(manhattan_dist_double);

    // This seems to not get compiled correctly? At least, it is cited as a
    // missing symbol at runtime
    template <typename Type>
    HWY_DLLEXPORT Type simd_manhattan_dist(
        const Type* x,
        const Type* y,
        const size_t size
    ){
        if(std::is_same<Type, float>::value){
            return HWY_DYNAMIC_DISPATCH(manhattan_dist_float)(x,  y, size);
        }
        else{
            return HWY_DYNAMIC_DISPATCH(manhattan_dist_double)(x,  y, size);
        }
    }
}
#endif  // HWY_ONCE
