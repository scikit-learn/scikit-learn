#include "simd.hpp"
/*
Externed declarations to provide separate translation unit for compilation
with specific feature flags (in this case -mavx)
*/
template float xsimd_manhattan_dist(const float* a, const float* b, const std::size_t size);
template double xsimd_manhattan_dist(const double* a, const double* b, const std::size_t size);
