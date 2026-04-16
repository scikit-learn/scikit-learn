cdef extern from "<complex.h>" nogil:
    # Trigonometric functions.
    double complex cacos(double complex z)
    double complex casin(double complex z)
    double complex catan(double complex z)
    double complex ccos(double complex z)
    double complex csin(double complex z)
    double complex ctan(double complex z)

    # Hyperbolic functions.
    double complex cacosh(double complex z)
    double complex casinh(double complex z)
    double complex catanh(double complex z)
    double complex ccosh(double complex z)
    double complex csinh(double complex z)
    double complex ctanh(double complex z)

    # Exponential and logarithmic functions.
    double complex cexp(double complex z)
    double complex clog(double complex z)
    double complex clog10(double complex z)

    # Power functions.
    double complex cpow(double complex x, double complex y)
    double complex csqrt(double complex z)

    # Absolute value, conjugates, and projection.
    double cabs(double complex z)
    double carg(double complex z)
    double complex conj(double complex z)
    double complex cproj(double complex z)

    # Decomposing complex values.
    double cimag(double complex z)
    double creal(double complex z)
