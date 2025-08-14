# NumPy math library
#
# This exports the functionality of the NumPy core math library, aka npymath,
# which provides implementations of C99 math functions and macros for system
# with a C89 library (such as MSVC). npymath is available with NumPy >=1.3,
# although some functions will require later versions. The spacing function is
# not in C99, but comes from Fortran.
#
# On the Cython side, the npymath functions are available without the "npy_"
# prefix that they have in C, to make this is a drop-in replacement for
# libc.math. The same is true for the constants, where possible.
#
# See the NumPy documentation for linking instructions.
#
# Complex number support and NumPy 2.0 half-precision functions are currently
# not exported.
#
# Author: Lars Buitinck


# NOTE: Deprecated since Cython 3.1.
#
# This file is mostly redundant with "libc.math" which should be used instead.


from libc.math cimport (
    NAN,
    INFINITY,

    isinf,
    isfinite,
    isnan,
    signbit,

    M_E as E,
    M_LOG2E as LOG2E,
    M_LOG10E as LOG10E,
    M_PI as PI,
    M_PI_2 as PI_2,
    M_PI_4 as PI_4,
    M_1_PI as NPY_1_PI,
    M_2_PI as NPY_2_PI,
    M_LN2 as LOGE2,
    M_LN10 as LOGE10,

    copysignf,
    nextafterf,
    copysign,
    nextafter,
    copysignl,
    nextafterl,

    sinf,
    cosf,
    tanf,
    sinhf,
    coshf,
    tanhf,
    fabsf,
    floorf,
    ceilf,
    rintf,
    sqrtf,
    log10f,
    logf,
    expf,
    expm1f,
    asinf,
    acosf,
    atanf,
    asinhf,
    acoshf,
    atanhf,
    log1pf,
    exp2f,
    log2f,
    atan2f,
    hypotf,
    powf,
    fmodf,
    modff,

    sinl,
    cosl,
    tanl,
    sinhl,
    coshl,
    tanhl,
    fabsl,
    floorl,
    ceill,
    rintl,
    sqrtl,
    log10l,
    logl,
    expl,
    expm1l,
    asinl,
    acosl,
    atanl,
    asinhl,
    acoshl,
    atanhl,
    log1pl,
    exp2l,
    log2l,
    atan2l,
    hypotl,
    powl,
    fmodl,
    modfl,
)


cdef extern from "numpy/npy_math.h" nogil:
    """
    #ifdef _MSC_VER
    #pragma message ("The 'numpy.math' import is outdated and deprecated. Use the standard 'libc.math' instead.")
    #else
    #warning The 'numpy.math' import is outdated and deprecated. Use the standard 'libc.math' instead.
    #endif
    """

    # Floating-point classification
    long double PZERO "NPY_PZERO"        # positive zero
    long double NZERO "NPY_NZERO"        # negative zero

    # Math constants
    long double EULER "NPY_EULER"       # Euler constant (gamma, 0.57721)

    # Low-level floating point manipulation (NumPy >=1.4)
    float spacingf "npy_spacingf"(float x)
    double spacing "npy_spacing"(double x)
    long double spacingl "npy_spacingl"(long double x)

    # NumPy extensions
    float deg2radf "npy_deg2radf"(float x)
    float rad2degf "npy_rad2degf"(float x)
    float logaddexpf "npy_logaddexpf"(float x, float y)
    float logaddexp2f "npy_logaddexp2f"(float x, float y)

    double deg2rad "npy_deg2rad"(double x)
    double rad2deg "npy_rad2deg"(double x)
    double logaddexp "npy_logaddexp"(double x, double y)
    double logaddexp2 "npy_logaddexp2"(double x, double y)

    long double deg2radl "npy_deg2radl"(long double x)
    long double rad2degl "npy_rad2degl"(long double x)
    long double logaddexpl "npy_logaddexpl"(long double x, long double y)
    long double logaddexp2l "npy_logaddexp2l"(long double x, long double y)
