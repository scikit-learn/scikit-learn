#include "xsimd.hpp"

namespace xs = xsimd;

/*Manhattan Distance*/
template <typename Type>
Type xsimd_manhattan_dist(const Type* a, const Type* b, const std::size_t size){
    using batch_type = xs::batch<Type, xs::avx>;

    // Unroll explicitly for pipeline optimization
    batch_type simd_x_0;
    batch_type simd_y_0;
    batch_type sum_0 = batch_type::broadcast((Type) 0.);

    batch_type simd_x_1;
    batch_type simd_y_1;
    batch_type sum_1 = batch_type::broadcast((Type) 0.);
    // End unrolled

    // Begin VECTOR LOOP
    std::size_t inc = batch_type::size;
    std::size_t loop_iter = inc * 2;
    std::size_t vec_size = size - size % loop_iter;
    std::size_t vec_remainder_size = size - size % inc;
    for(std::size_t idx = 0; idx < vec_size; idx += loop_iter) {
        // Begin unrolled
        simd_x_0 = batch_type::load_unaligned(&a[idx + inc * 0]);
        simd_y_0 = batch_type::load_unaligned(&b[idx + inc * 0]);
        sum_0 += xs::fabs(simd_x_0 - simd_y_0);

        simd_x_1 = batch_type::load_unaligned(&a[idx + inc * 1]);
        simd_y_1 = batch_type::load_unaligned(&b[idx + inc * 1]);
        sum_1 += xs::fabs(simd_x_1 - simd_y_1);
        // End unrolled

    }
    for(std::size_t idx = vec_size; idx < vec_remainder_size; idx += inc) {
        simd_x_0 = batch_type::load_unaligned(&a[idx + inc * 0]);
        simd_y_0 = batch_type::load_unaligned(&b[idx + inc * 0]);
        sum_0 += xs::fabs(simd_x_0 - simd_y_0);
    }
    // End VECTOR LOOP
    // Reduction
    sum_0 += sum_1;
    batch_type batch_sum = xs::reduce_add(sum_0);
    double scalar_sum = *(Type*)&batch_sum;

    // Remaining scalar loop that cannot be vectorized
    for(std::size_t idx = vec_remainder_size; idx < size; ++idx) {
        scalar_sum += fabs(a[idx] - b[idx]);
    }
    return (Type) scalar_sum;
}

/*Extern declarations to later be mapped to simd.cpp declarations at link time*/
extern template float xsimd_manhattan_dist(const float* a, const float* b, const std::size_t size);
extern template double xsimd_manhattan_dist(const double* a, const double* b, const std::size_t size);
