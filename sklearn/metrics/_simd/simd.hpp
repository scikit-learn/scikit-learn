#include <hwy/highway.h>
#include <iostream>
#include <cmath>
namespace hn = hwy::HWY_NAMESPACE;

template <typename Type>
Type simd_manhattan_dist(
    const Type* x,
    const Type* y,
    const size_t size
) {
    const hn::ScalableTag<Type> d;
    auto simd_sum_1 = hn::Zero(d);
    auto simd_sum_2 = hn::Zero(d);
    auto lane_step = hn::Lanes(d);
    size_t loop_iter = lane_step * 2;
    size_t vec_size = size - size % loop_iter;
    size_t vec_remainder_size = size - size % lane_step;

    for (size_t i = 0; i < vec_size; i += loop_iter) {
        const auto simd_x_1 = hn::Load(d, x + i);
        const auto simd_y_1 = hn::Load(d, y + i);
        simd_sum_1 += hn::AbsDiff(simd_x_1, simd_y_1);

        const auto simd_x_2 = hn::Load(d, x + i + lane_step);
        const auto simd_y_2 = hn::Load(d, y + i + lane_step);
        simd_sum_2 += hn::AbsDiff(simd_x_2, simd_y_2);
    }
    for (size_t i = vec_size; i < vec_remainder_size; i += lane_step) {
        const auto simd_x_1 = hn::Load(d, x + i);
        const auto simd_y_1 = hn::Load(d, y + i);
        simd_sum_1 += hn::AbsDiff(simd_x_1, simd_y_1);
    }
    simd_sum_1 += simd_sum_2;
    Type scalar_sum = hn::ReduceSum(d, simd_sum_1);
    for (size_t i = vec_remainder_size; i < size; i += 1) {
        scalar_sum += fabs(x[i] - y[i]);
    }
    return scalar_sum;
}

/*Extern declarations to later be mapped to simd.cpp declarations at link time*/
extern template float simd_manhattan_dist(
    const float* x,
    const float* y,
    const size_t size
);
extern template double simd_manhattan_dist(
    const double* x,
    const double* y,
    const size_t size
);
