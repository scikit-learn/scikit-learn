
cdef extern from "<cmath>" namespace "std" nogil:
    # all C99 functions
    float acos(float x) except +
    double acos(double x) except +
    long double acos(long double x) except +
    float acosf(float x) except +
    long double acosl(long double x) except +

    float asin(float x) except +
    double asin(double x) except +
    long double asin(long double x) except +
    float asinf(float x) except +
    long double asinl(long double x) except +

    float atan(float x) except +
    double atan(double x) except +
    long double atan(long double x) except +
    float atanf(float x) except +
    long double atanl(long double x) except +

    float atan2(float y, float x) except +
    double atan2(double y, double x) except +
    long double atan2(long double y, long double x) except +
    float atan2f(float y, float x) except +
    long double atan2l(long double y, long double x) except +

    float cos(float x) except +
    double cos(double x) except +
    long double cos(long double x) except +
    float cosf(float x) except +
    long double cosl(long double x) except +

    float sin(float x) except +
    double sin(double x) except +
    long double sin(long double x) except +
    float sinf(float x) except +
    long double sinl(long double x) except +

    float tan(float x) except +
    double tan(double x) except +
    long double tan(long double x) except +
    float tanf(float x) except +
    long double tanl(long double x) except +

    float acosh(float x) except +
    double acosh(double x) except +
    long double acosh(long double x) except +
    float acoshf(float x) except +
    long double acoshl(long double x) except +

    float asinh(float x) except +
    double asinh(double x) except +
    long double asinh(long double x) except +
    float asinhf(float x) except +
    long double asinhl(long double x) except +

    float atanh(float x) except +
    double atanh(double x) except +
    long double atanh(long double x) except +
    float atanhf(float x) except +
    long double atanhl(long double x) except +

    float cosh(float x) except +
    double cosh(double x) except +
    long double cosh(long double x) except +
    float coshf(float x) except +
    long double coshl(long double x) except +

    float sinh(float x) except +
    double sinh(double x) except +
    long double sinh(long double x) except +
    float sinhf(float x) except +
    long double sinhl(long double x) except +

    float tanh(float x) except +
    double tanh(double x) except +
    long double tanh(long double x) except +
    float tanhf(float x) except +
    long double tanhl(long double x) except +

    float exp(float x) except +
    double exp(double x) except +
    long double exp(long double x) except +
    float expf(float x) except +
    long double expl(long double x) except +

    float exp2(float x) except +
    double exp2(double x) except +
    long double exp2(long double x) except +
    float exp2f(float x) except +
    long double exp2l(long double x) except +

    float expm1(float x) except +
    double expm1(double x) except +
    long double expm1(long double x) except +
    float expm1f(float x) except +
    long double expm1l(long double x) except +

    float frexp(float value, int* exp) except +
    double frexp(double value, int* exp) except +
    long double frexp(long double value, int* exp) except +
    float frexpf(float value, int* exp) except +
    long double frexpl(long double value, int* exp) except +

    int ilogb(float x) except +
    int ilogb(double x) except +
    int ilogb(long double x) except +
    int ilogbf(float x) except +
    int ilogbl(long double x) except +

    float ldexp(float x, int exp) except +
    double ldexp(double x, int exp) except +
    long double ldexp(long double x, int exp) except +
    float ldexpf(float x, int exp) except +
    long double ldexpl(long double x, int exp) except +

    float log(float x) except +
    double log(double x) except +
    long double log(long double x) except +
    float logf(float x) except +
    long double logl(long double x) except +

    float log10(float x) except +
    double log10(double x) except +
    long double log10(long double x) except +
    float log10f(float x) except +
    long double log10l(long double x) except +

    float log1p(float x) except +
    double log1p(double x) except +
    long double log1p(long double x) except +
    float log1pf(float x) except +
    long double log1pl(long double x) except +

    float log2(float x) except +
    double log2(double x) except +
    long double log2(long double x) except +
    float log2f(float x) except +
    long double log2l(long double x) except +

    float logb(float x) except +
    double logb(double x) except +
    long double logb(long double x) except +
    float logbf(float x) except +
    long double logbl(long double x) except +

    float modf(float value, float* iptr) except +
    double modf(double value, double* iptr) except +
    long double modf(long double value, long double* iptr) except +
    float modff(float value, float* iptr) except +
    long double modfl(long double value, long double* iptr) except +

    float scalbn(float x, int n) except +
    double scalbn(double x, int n) except +
    long double scalbn(long double x, int n) except +
    float scalbnf(float x, int n) except +
    long double scalbnl(long double x, int n) except +

    float scalbln(float x, long int n) except +
    double scalbln(double x, long int n) except +
    long double scalbln(long double x, long int n) except +
    float scalblnf(float x, long int n) except +
    long double scalblnl(long double x, long int n) except +

    float cbrt(float x) except +
    double cbrt(double x) except +
    long double cbrt(long double x) except +
    float cbrtf(float x) except +
    long double cbrtl(long double x) except +

    # absolute values
    int abs(int j) except +
    long int abs(long int j) except +
    long long int abs(long long int j) except +
    float abs(float j) except +
    double abs(double j) except +
    long double abs(long double j) except +

    float fabs(float x) except +
    double fabs(double x) except +
    long double fabs(long double x) except +
    float fabsf(float x) except +
    long double fabsl(long double x) except +

    float hypot(float x, float y) except +
    double hypot(double x, double y) except +
    long double hypot(long double x, long double y) except +
    float hypotf(float x, float y) except +
    long double hypotl(long double x, long double y) except +

    # C++17 three-dimensional hypotenuse
    float hypot(float x, float y, float z) except +
    double hypot(double x, double y, double z) except +
    long double hypot(long double x, long double y, long double z) except +

    float pow(float x, float y) except +
    double pow(double x, double y) except +
    long double pow(long double x, long double y) except +
    float powf(float x, float y) except +
    long double powl(long double x, long double y) except +

    float sqrt(float x) except +
    double sqrt(double x) except +
    long double sqrt(long double x) except +
    float sqrtf(float x) except +
    long double sqrtl(long double x) except +

    float erf(float x) except +
    double erf(double x) except +
    long double erf(long double x) except +
    float erff(float x) except +
    long double erfl(long double x) except +

    float erfc(float x) except +
    double erfc(double x) except +
    long double erfc(long double x) except +
    float erfcf(float x) except +
    long double erfcl(long double x) except +

    float lgamma(float x) except +
    double lgamma(double x) except +
    long double lgamma(long double x) except +
    float lgammaf(float x) except +
    long double lgammal(long double x) except +

    float tgamma(float x) except +
    double tgamma(double x) except +
    long double tgamma(long double x) except +
    float tgammaf(float x) except +
    long double tgammal(long double x) except +

    float ceil(float x) except +
    double ceil(double x) except +
    long double ceil(long double x) except +
    float ceilf(float x) except +
    long double ceill(long double x) except +

    float floor(float x) except +
    double floor(double x) except +
    long double floor(long double x) except +
    float floorf(float x) except +
    long double floorl(long double x) except +

    float nearbyint(float x) except +
    double nearbyint(double x) except +
    long double nearbyint(long double x) except +
    float nearbyintf(float x) except +
    long double nearbyintl(long double x) except +

    float rint(float x) except +
    double rint(double x) except +
    long double rint(long double x) except +
    float rintf(float x) except +
    long double rintl(long double x) except +

    long int lrint(float x) except +
    long int lrint(double x) except +
    long int lrint(long double x) except +
    long int lrintf(float x) except +
    long int lrintl(long double x) except +

    long long int llrint(float x) except +
    long long int llrint(double x) except +
    long long int llrint(long double x) except +
    long long int llrintf(float x) except +
    long long int llrintl(long double x) except +

    float round(float x) except +
    double round(double x) except +
    long double round(long double x) except +
    float roundf(float x) except +
    long double roundl(long double x) except +

    long int lround(float x) except +
    long int lround(double x) except +
    long int lround(long double x) except +
    long int lroundf(float x) except +
    long int lroundl(long double x) except +

    long long int llround(float x) except +
    long long int llround(double x) except +
    long long int llround(long double x) except +
    long long int llroundf(float x) except +
    long long int llroundl(long double x) except +

    float trunc(float x) except +
    double trunc(double x) except +
    long double trunc(long double x) except +
    float truncf(float x) except +
    long double truncl(long double x) except +

    float fmod(float x, float y) except +
    double fmod(double x, double y) except +
    long double fmod(long double x, long double y) except +
    float fmodf(float x, float y) except +
    long double fmodl(long double x, long double y) except +

    float remainder(float x, float y) except +
    double remainder(double x, double y) except +
    long double remainder(long double x, long double y) except +
    float remainderf(float x, float y) except +
    long double remainderl(long double x, long double y) except +

    float remquo(float x, float y, int* quo) except +
    double remquo(double x, double y, int* quo) except +
    long double remquo(long double x, long double y, int* quo) except +
    float remquof(float x, float y, int* quo) except +
    long double remquol(long double x, long double y, int* quo) except +

    float copysign(float x, float y) except +
    double copysign(double x, double y) except +
    long double copysign(long double x, long double y) except +
    float copysignf(float x, float y) except +
    long double copysignl(long double x, long double y) except +

    double nan(const char* tagp) except +
    float nanf(const char* tagp) except +
    long double nanl(const char* tagp) except +

    float nextafter(float x, float y) except +
    double nextafter(double x, double y) except +
    long double nextafter(long double x, long double y) except +
    float nextafterf(float x, float y) except +
    long double nextafterl(long double x, long double y) except +

    float nexttoward(float x, long double y) except +
    double nexttoward(double x, long double y) except +
    long double nexttoward(long double x, long double y) except +
    float nexttowardf(float x, long double y) except +
    long double nexttowardl(long double x, long double y) except +

    float fdim(float x, float y) except +
    double fdim(double x, double y) except +
    long double fdim(long double x, long double y) except +
    float fdimf(float x, float y) except +
    long double fdiml(long double x, long double y) except +

    float fmax(float x, float y) except +
    double fmax(double x, double y) except +
    long double fmax(long double x, long double y) except +
    float fmaxf(float x, float y) except +
    long double fmaxl(long double x, long double y) except +

    float fmin(float x, float y) except +
    double fmin(double x, double y) except +
    long double fmin(long double x, long double y) except +
    float fminf(float x, float y) except +
    long double fminl(long double x, long double y) except +

    float fma(float x, float y, float z) except +
    double fma(double x, double y, double z) except +
    long double fma(long double x, long double y, long double z) except +
    float fmaf(float x, float y, float z) except +
    long double fmal(long double x, long double y, long double z) except +

    # C++20 linear interpolation
    float lerp(float a, float b, float t)
    double lerp(double a, double b, double t)
    long double lerp(long double a, long double b, long double t)

    # classification / comparison functions
    int fpclassify(float x) except +
    int fpclassify(double x) except +
    int fpclassify(long double x) except +

    bint isfinite(float x) except +
    bint isfinite(double x) except +
    bint isfinite(long double x) except +

    bint isinf(float x) except +
    bint isinf(double x) except +
    bint isinf(long double x) except +

    bint isnan(float x) except +
    bint isnan(double x) except +
    bint isnan(long double x) except +

    bint isnormal(float x) except +
    bint isnormal(double x) except +
    bint isnormal(long double x) except +

    bint signbit(float x) except +
    bint signbit(double x) except +
    bint signbit(long double x) except +

    bint isgreater(float x, float y) except +
    bint isgreater(double x, double y) except +
    bint isgreater(long double x, long double y) except +

    bint isgreaterequal(float x, float y) except +
    bint isgreaterequal(double x, double y) except +
    bint isgreaterequal(long double x, long double y) except +

    bint isless(float x, float y) except +
    bint isless(double x, double y) except +
    bint isless(long double x, long double y) except +

    bint islessequal(float x, float y) except +
    bint islessequal(double x, double y) except +
    bint islessequal(long double x, long double y) except +

    bint islessgreater(float x, float y) except +
    bint islessgreater(double x, double y) except +
    bint islessgreater(long double x, long double y) except +

    bint isunordered(float x, float y) except +
    bint isunordered(double x, double y) except +
    bint isunordered(long double x, long double y) except +

    # C++17 mathematical special functions

    # associated Laguerre polynomials
    double       assoc_laguerre(unsigned int n, unsigned int m, double x) except +
    float        assoc_laguerref(unsigned int n, unsigned int m, float x) except +
    long double  assoc_laguerrel(unsigned int n, unsigned int m, long double x) except +

    # associated Legendre functions
    double       assoc_legendre(unsigned int l, unsigned int m, double x) except +
    float        assoc_legendref(unsigned int l, unsigned int m, float x) except +
    long double  assoc_legendrel(unsigned int l, unsigned int m, long double x) except +

    # beta function
    double       beta(double x, double y) except +
    float        betaf(float x, float y) except +
    long double  betal(long double x, long double y) except +

    # complete elliptic integral of the first kind
    double       comp_ellint_1(double k) except +
    float        comp_ellint_1f(float k) except +
    long double  comp_ellint_1l(long double k) except +

    # complete elliptic integral of the second kind
    double       comp_ellint_2(double k) except +
    float        comp_ellint_2f(float k) except +
    long double  comp_ellint_2l(long double k) except +

    # complete elliptic integral of the third kind
    double       comp_ellint_3(double k, double nu) except +
    float        comp_ellint_3f(float k, float nu) except +
    long double  comp_ellint_3l(long double k, long double nu) except +

    # regular modified cylindrical Bessel functions
    double       cyl_bessel_i(double nu, double x) except +
    float        cyl_bessel_if(float nu, float x) except +
    long double  cyl_bessel_il(long double nu, long double x) except +

    # cylindrical Bessel functions of the first kind
    double       cyl_bessel_j(double nu, double x) except +
    float        cyl_bessel_jf(float nu, float x) except +
    long double  cyl_bessel_jl(long double nu, long double x) except +

    # irregular modified cylindrical Bessel functions
    double       cyl_bessel_k(double nu, double x) except +
    float        cyl_bessel_kf(float nu, float x) except +
    long double  cyl_bessel_kl(long double nu, long double x) except +

    # cylindrical Neumann functions
    # cylindrical Bessel functions of the second kind
    double       cyl_neumann(double nu, double x) except +
    float        cyl_neumannf(float nu, float x) except +
    long double  cyl_neumannl(long double nu, long double x) except +

    # incomplete elliptic integral of the first kind
    double       ellint_1(double k, double phi) except +
    float        ellint_1f(float k, float phi) except +
    long double  ellint_1l(long double k, long double phi) except +

    # incomplete elliptic integral of the second kind
    double       ellint_2(double k, double phi) except +
    float        ellint_2f(float k, float phi) except +
    long double  ellint_2l(long double k, long double phi) except +

    # incomplete elliptic integral of the third kind
    double       ellint_3(double k, double nu, double phi) except +
    float        ellint_3f(float k, float nu, float phi) except +
    long double  ellint_3l(long double k, long double nu, long double phi) except +

    # exponential integral
    double       expint(double x) except +
    float        expintf(float x) except +
    long double  expintl(long double x) except +

    # Hermite polynomials
    double       hermite(unsigned int n, double x) except +
    float        hermitef(unsigned int n, float x) except +
    long double  hermitel(unsigned int n, long double x) except +

    # Laguerre polynomials
    double       laguerre(unsigned int n, double x) except +
    float        laguerref(unsigned int n, float x) except +
    long double  laguerrel(unsigned int n, long double x) except +

    # Legendre polynomials
    double       legendre(unsigned int l, double x) except +
    float        legendref(unsigned int l, float x) except +
    long double  legendrel(unsigned int l, long double x) except +

    # Riemann zeta function
    double       riemann_zeta(double x) except +
    float        riemann_zetaf(float x) except +
    long double  riemann_zetal(long double x) except +

    # spherical Bessel functions of the first kind
    double       sph_bessel(unsigned int n, double x) except +
    float        sph_besself(unsigned int n, float x) except +
    long double  sph_bessell(unsigned int n, long double x) except +

    # spherical associated Legendre functions
    double       sph_legendre(unsigned int l, unsigned int m, double theta) except +
    float        sph_legendref(unsigned int l, unsigned int m, float theta) except +
    long double  sph_legendrel(unsigned int l, unsigned int m, long double theta) except +

    # spherical Neumann functions
    # spherical Bessel functions of the second kind
    double       sph_neumann(unsigned int n, double x) except +
    float        sph_neumannf(unsigned int n, float x) except +
    long double  sph_neumannl(unsigned int n, long double x) except +
