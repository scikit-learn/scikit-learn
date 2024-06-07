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
#define SPECFUN_HOST_DEVICE __host__ __device__

#include <cuda/std/cmath>
#include <cuda/std/limits>

// Fallback to global namespace for functions unsupported on NVRTC Jit
#ifdef _LIBCUDACXX_COMPILER_NVRTC
#include <cuda_runtime.h>
#endif

namespace std {

SPECFUN_HOST_DEVICE inline double abs(double num) { return cuda::std::abs(num); }

SPECFUN_HOST_DEVICE inline double exp(double num) { return cuda::std::exp(num); }

SPECFUN_HOST_DEVICE inline double log(double num) { return cuda::std::log(num); }

SPECFUN_HOST_DEVICE inline double sqrt(double num) { return cuda::std::sqrt(num); }

SPECFUN_HOST_DEVICE inline bool isnan(double num) { return cuda::std::isnan(num); }

SPECFUN_HOST_DEVICE inline bool isfinite(double num) { return cuda::std::isfinite(num); }

SPECFUN_HOST_DEVICE inline double pow(double x, double y) { return cuda::std::pow(x, y); }

SPECFUN_HOST_DEVICE inline double sin(double x) { return cuda::std::sin(x); }

SPECFUN_HOST_DEVICE inline double tan(double x) { return cuda::std::tan(x); }

SPECFUN_HOST_DEVICE inline double sinh(double x) { return cuda::std::sinh(x); }

SPECFUN_HOST_DEVICE inline double cosh(double x) { return cuda::std::cosh(x); }

SPECFUN_HOST_DEVICE inline bool signbit(double x) { return cuda::std::signbit(x); }

// Fallback to global namespace for functions unsupported on NVRTC
#ifndef _LIBCUDACXX_COMPILER_NVRTC
SPECFUN_HOST_DEVICE inline double ceil(double x) { return cuda::std::ceil(x); }
SPECFUN_HOST_DEVICE inline double floor(double x) { return cuda::std::floor(x); }
SPECFUN_HOST_DEVICE inline double trunc(double x) { return cuda::std::trunc(x); }
SPECFUN_HOST_DEVICE inline double fma(double x, double y, double z) { return cuda::std::fma(x, y, z); }
SPECFUN_HOST_DEVICE inline double copysign(double x, double y) { return cuda::std::copysign(x, y); }
SPECFUN_HOST_DEVICE inline double modf(double value, double *iptr) { return cuda::std::modf(value, iptr); }

#else
SPECFUN_HOST_DEVICE inline double ceil(double x) { return ::ceil(x); }
SPECFUN_HOST_DEVICE inline double floor(double x) { return ::floor(x); }
SPECFUN_HOST_DEVICE inline double trunc(double x) { return ::trunc(x); }
SPECFUN_HOST_DEVICE inline double fma(double x, double y, double z) { return ::fma(x, y, z); }
SPECFUN_HOST_DEVICE inline double copysign(double x, double y) { return ::copysign(x, y); }
SPECFUN_HOST_DEVICE inline double modf(double value, double *iptr) { return ::modf(value, iptr); }
#endif

template <typename T>
using numeric_limits = cuda::std::numeric_limits<T>;

// Must use thrust for complex types in order to support CuPy
template <typename T>
using complex = thrust::complex<T>;

template <typename T>
SPECFUN_HOST_DEVICE T abs(const complex<T> &z) {
    return thrust::abs(z);
}

template <typename T>
SPECFUN_HOST_DEVICE complex<T> exp(const complex<T> &z) {
    return thrust::exp(z);
}

template <typename T>
SPECFUN_HOST_DEVICE complex<T> log(const complex<T> &z) {
    return thrust::log(z);
}

template <typename T>
SPECFUN_HOST_DEVICE T norm(const complex<T> &z) {
    return thrust::norm(z);
}

template <typename T>
SPECFUN_HOST_DEVICE complex<T> sqrt(const complex<T> &z) {
    return thrust::sqrt(z);
}

template <typename T>
SPECFUN_HOST_DEVICE complex<T> conj(const complex<T> &z) {
    return thrust::conj(z);
}

} // namespace std

#else
#define SPECFUN_HOST_DEVICE

#include <cmath>
#include <complex>
#include <limits>
#include <math.h>

#endif
