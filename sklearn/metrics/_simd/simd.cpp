// #define HWY_TARGETS HWY_AVX2
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
            const hn::ScalableTag<Type> d;
            using batch_type = decltype(hn::Zero(d));

            batch_type simd_sum_1 = hn::Zero(d);
            batch_type simd_x_1;
            batch_type simd_y_1;

            batch_type simd_sum_2 = hn::Zero(d);
            batch_type simd_x_2;
            batch_type simd_y_2;

            size_t lane_size = hn::Lanes(d);
            size_t loop_iter = lane_size * 2;
            size_t vec_size = size - size % loop_iter;
            size_t vec_remainder_size = size - size % lane_size;

            for (size_t i = 0; i < vec_size; i += loop_iter) {
                simd_x_1 = hn::LoadU(d, x + i);
                simd_y_1 = hn::LoadU(d, y + i);
                simd_sum_1 += hn::AbsDiff(simd_x_1, simd_y_1);

                simd_x_2 = hn::LoadU(d, x + i + lane_size);
                simd_y_2 = hn::LoadU(d, y + i + lane_size);
                simd_sum_2 += hn::AbsDiff(simd_x_2, simd_y_2);
            }
            for (size_t i = vec_size; i < vec_remainder_size; i += loop_iter) {
                simd_x_1 = hn::LoadU(d, x + i);
                simd_y_1 = hn::LoadU(d, y + i);
                simd_sum_1 += hn::AbsDiff(simd_x_1, simd_y_1);
            }
            simd_sum_1 += simd_sum_2;
            Type scalar_sum = hn::ReduceSum(d, simd_sum_1);
            for (size_t i = vec_remainder_size; i < size; i += 1) {
                scalar_sum += fabs(x[i] - y[i]);
            }
            return scalar_sum;
        }
        float manhattan_dist_f32(
            const float* x,
            const float* y,
            const size_t size
        ) {
            return manhattan_dist<float>(x, y, size);
        }
        double manhattan_dist_f64(
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
#include "simd.hpp"

namespace manhattan {

    HWY_EXPORT(manhattan_dist_f32);
    HWY_EXPORT(manhattan_dist_f64);

    HWY_DLLEXPORT float simd_manhattan_dist_f32(
        const float* x,
        const float* y,
        const size_t size
    ){
        return HWY_DYNAMIC_DISPATCH(manhattan_dist_f32)(x,  y, size);
    }
    HWY_DLLEXPORT double simd_manhattan_dist_f64(
        const double* x,
        const double* y,
        const size_t size
    ){
        return HWY_DYNAMIC_DISPATCH(manhattan_dist_f64)(x,  y, size);
    }
}
#endif  // HWY_ONCE
