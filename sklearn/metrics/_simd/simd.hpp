#include <highway.h>
namespace hn = hwy::HWY_NAMESPACE;

template <typename Type>
Type simd_manhattan_dist(
    const Type* HWY_RESTRICT x,
    const Type* HWY_RESTRICT y,
    const size_t size,
) {
  const hn::ScalableTag<Type> d;
  auto simd_sum_1 = hn::Zero(d);
  auto simd_sum_2 = hn::Zero(d);
  auto lane_step = hn::Lanes(d);
  for (size_t i = 0; i < size; i += 2 * lane_step) {
    const auto simd_x_1 = hn::Load(d, x + i);
    const auto simd_y_1 = hn::Load(d, y + i);
    simd_sum_1 += hn::AbsDiff(simd_x_1, simd_y_1);

    const auto simd_x_2 = hn::Load(d, x + i + lane_step);
    const auto simd_y_2 = hn::Load(d, y + i + lane_step);
    simd_sum_2 += hn::AbsDiff(simd_x_2, simd_y_2);
  }
  simd_sum_1 += simd_sum_2;
  return hn::ReduceSum(simd_sum_1);
}

/*Extern declarations to later be mapped to simd.cpp declarations at link time*/
extern template float simd_manhattan_dist(
    const float* HWY_RESTRICT x,
    const float* HWY_RESTRICT y,
    const size_t size
);
extern template double simd_manhattan_dist(
    const double* HWY_RESTRICT x,
    const double* HWY_RESTRICT y,
    const size_t size
);
