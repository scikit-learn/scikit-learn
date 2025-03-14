#pragma once

// Define math constants if they are not available
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#ifndef M_LOG2E
#define M_LOG2E 1.44269504088896340736
#endif

#ifndef M_LOG10E
#define M_LOG10E 0.434294481903251827651
#endif

#ifndef M_LN2
#define M_LN2 0.693147180559945309417
#endif

#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI_4
#define M_PI_4 0.785398163397448309616
#endif

#ifndef M_1_PI
#define M_1_PI 0.318309886183790671538
#endif

#ifndef M_2_PI
#define M_2_PI 0.636619772367581343076
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524401
#endif

#ifdef __CUDACC__
#define XSF_HOST_DEVICE __host__ __device__

#include <cuda/std/cmath>
#include <cuda/std/cstddef>
#include <cuda/std/cstdint>
#include <cuda/std/limits>
#include <cuda/std/tuple>
#include <cuda/std/type_traits>
#include <cuda/std/utility>

// Fallback to global namespace for functions unsupported on NVRTC Jit
#ifdef _LIBCUDACXX_COMPILER_NVRTC
#include <cuda_runtime.h>
#endif

namespace std {

XSF_HOST_DEVICE inline double abs(double num) { return cuda::std::abs(num); }

XSF_HOST_DEVICE inline double exp(double num) { return cuda::std::exp(num); }

XSF_HOST_DEVICE inline double log(double num) { return cuda::std::log(num); }

XSF_HOST_DEVICE inline double sqrt(double num) { return cuda::std::sqrt(num); }

XSF_HOST_DEVICE inline bool isinf(double num) { return cuda::std::isinf(num); }

XSF_HOST_DEVICE inline bool isnan(double num) { return cuda::std::isnan(num); }

XSF_HOST_DEVICE inline bool isfinite(double num) { return cuda::std::isfinite(num); }

XSF_HOST_DEVICE inline double pow(double x, double y) { return cuda::std::pow(x, y); }

XSF_HOST_DEVICE inline double sin(double x) { return cuda::std::sin(x); }

XSF_HOST_DEVICE inline double cos(double x) { return cuda::std::cos(x); }

XSF_HOST_DEVICE inline double tan(double x) { return cuda::std::tan(x); }

XSF_HOST_DEVICE inline double atan(double x) { return cuda::std::atan(x); }

XSF_HOST_DEVICE inline double acos(double x) { return cuda::std::acos(x); }

XSF_HOST_DEVICE inline double sinh(double x) { return cuda::std::sinh(x); }

XSF_HOST_DEVICE inline double cosh(double x) { return cuda::std::cosh(x); }

XSF_HOST_DEVICE inline double asinh(double x) { return cuda::std::asinh(x); }

XSF_HOST_DEVICE inline bool signbit(double x) { return cuda::std::signbit(x); }

// Fallback to global namespace for functions unsupported on NVRTC
#ifndef _LIBCUDACXX_COMPILER_NVRTC
XSF_HOST_DEVICE inline double ceil(double x) { return cuda::std::ceil(x); }
XSF_HOST_DEVICE inline double floor(double x) { return cuda::std::floor(x); }
XSF_HOST_DEVICE inline double round(double x) { return cuda::std::round(x); }
XSF_HOST_DEVICE inline double trunc(double x) { return cuda::std::trunc(x); }
XSF_HOST_DEVICE inline double fma(double x, double y, double z) { return cuda::std::fma(x, y, z); }
XSF_HOST_DEVICE inline double copysign(double x, double y) { return cuda::std::copysign(x, y); }
XSF_HOST_DEVICE inline double modf(double value, double *iptr) { return cuda::std::modf(value, iptr); }
XSF_HOST_DEVICE inline double fmax(double x, double y) { return cuda::std::fmax(x, y); }
XSF_HOST_DEVICE inline double fmin(double x, double y) { return cuda::std::fmin(x, y); }
XSF_HOST_DEVICE inline double log10(double num) { return cuda::std::log10(num); }
XSF_HOST_DEVICE inline double log1p(double num) { return cuda::std::log1p(num); }
XSF_HOST_DEVICE inline double frexp(double num, int *exp) { return cuda::std::frexp(num, exp); }
XSF_HOST_DEVICE inline double ldexp(double num, int exp) { return cuda::std::ldexp(num, exp); }
XSF_HOST_DEVICE inline double fmod(double x, double y) { return cuda::std::fmod(x, y); }
XSF_HOST_DEVICE inline double nextafter(double from, double to) { return cuda::std::nextafter(from, to); }
#else
XSF_HOST_DEVICE inline double ceil(double x) { return ::ceil(x); }
XSF_HOST_DEVICE inline double floor(double x) { return ::floor(x); }
XSF_HOST_DEVICE inline double round(double x) { return ::round(x); }
XSF_HOST_DEVICE inline double trunc(double x) { return ::trunc(x); }
XSF_HOST_DEVICE inline double fma(double x, double y, double z) { return ::fma(x, y, z); }
XSF_HOST_DEVICE inline double copysign(double x, double y) { return ::copysign(x, y); }
XSF_HOST_DEVICE inline double modf(double value, double *iptr) { return ::modf(value, iptr); }
XSF_HOST_DEVICE inline double fmax(double x, double y) { return ::fmax(x, y); }
XSF_HOST_DEVICE inline double fmin(double x, double y) { return ::fmin(x, y); }
XSF_HOST_DEVICE inline double log10(double num) { return ::log10(num); }
XSF_HOST_DEVICE inline double log1p(double num) { return ::log1p(num); }
XSF_HOST_DEVICE inline double frexp(double num, int *exp) { return ::frexp(num, exp); }
XSF_HOST_DEVICE inline double ldexp(double num, int exp) { return ::ldexp(num, exp); }
XSF_HOST_DEVICE inline double fmod(double x, double y) { return ::fmod(x, y); }
XSF_HOST_DEVICE inline double nextafter(double from, double to) { return ::nextafter(from, to); }
#endif

template <typename T>
XSF_HOST_DEVICE void swap(T &a, T &b) {
    cuda::std::swap(a, b);
}

// Reimplement std::clamp until it's available in CuPy
template <typename T>
XSF_HOST_DEVICE constexpr T clamp(T &v, T &lo, T &hi) {
    return v < lo ? lo : (v > hi ? lo : v);
}

template <typename T>
using numeric_limits = cuda::std::numeric_limits<T>;

// Must use thrust for complex types in order to support CuPy
template <typename T>
using complex = thrust::complex<T>;

template <typename T>
XSF_HOST_DEVICE T abs(const complex<T> &z) {
    return thrust::abs(z);
}

template <typename T>
XSF_HOST_DEVICE complex<T> exp(const complex<T> &z) {
    return thrust::exp(z);
}

template <typename T>
XSF_HOST_DEVICE complex<T> log(const complex<T> &z) {
    return thrust::log(z);
}

template <typename T>
XSF_HOST_DEVICE T norm(const complex<T> &z) {
    return thrust::norm(z);
}

template <typename T>
XSF_HOST_DEVICE complex<T> sqrt(const complex<T> &z) {
    return thrust::sqrt(z);
}

template <typename T>
XSF_HOST_DEVICE complex<T> conj(const complex<T> &z) {
    return thrust::conj(z);
}

template <typename T>
XSF_HOST_DEVICE complex<T> pow(const complex<T> &x, const complex<T> &y) {
    return thrust::pow(x, y);
}

template <typename T>
XSF_HOST_DEVICE complex<T> pow(const complex<T> &x, const T &y) {
    return thrust::pow(x, y);
}

// Other types and utilities
template <typename T>
using is_floating_point = cuda::std::is_floating_point<T>;

template <bool Cond, typename T = void>
using enable_if = cuda::std::enable_if<Cond, T>;

template <typename T>
using decay = cuda::std::decay<T>;

template <typename T>
using invoke_result = cuda::std::invoke_result<T>;

template <typename T1, typename T2>
using pair = cuda::std::pair<T1, T2>;

template <typename... Types>
using tuple = cuda::std::tuple<Types...>;

using cuda::std::ptrdiff_t;
using cuda::std::size_t;
using cuda::std::uint64_t;

#define XSF_ASSERT(a)

} // namespace std

#else
#define XSF_HOST_DEVICE

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <math.h>
#include <tuple>
#include <type_traits>
#include <utility>

#ifdef DEBUG
#define XSF_ASSERT(a) assert(a)
#else
#define XSF_ASSERT(a)
#endif

namespace xsf {

// basic
using std::abs;

// exponential
using std::exp;

// power
using std::sqrt;

// trigonometric
using std::cos;
using std::sin;

// floating-point manipulation
using std::copysign;

// classification and comparison
using std::isfinite;
using std::isinf;
using std::isnan;
using std::signbit;

// complex
using std::imag;
using std::real;

template <typename T>
struct remove_complex {
    using type = T;
};

template <typename T>
struct remove_complex<std::complex<T>> {
    using type = T;
};

template <typename T>
using remove_complex_t = typename remove_complex<T>::type;

template <typename T>
struct complex_type {
    using type = std::complex<T>;
};

template <typename T>
using complex_type_t = typename complex_type<T>::type;

template <typename T>
using complex = complex_type_t<T>;

} // namespace xsf

#endif
