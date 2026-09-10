// Thin extern "C" wrapper around x86-simd-sort's keyvalue_qsort template,
// instantiated for the (float32, int32) and (float64, int32) pairs used by
// sklearn.utils._sorting.simultaneous_sort. See _sorting.pyx.
#include "x86simdsort.h"

#include <cstddef>
#include <cstdint>

extern "C" {

void xss_keyvalue_sort_f64_i32(double *keys, int32_t *vals, size_t n) {
    x86simdsort::keyvalue_qsort<double, int32_t>(keys, vals, n, false, false);
}

void xss_keyvalue_sort_f32_i32(float *keys, int32_t *vals, size_t n) {
    x86simdsort::keyvalue_qsort<float, int32_t>(keys, vals, n, false, false);
}

}
