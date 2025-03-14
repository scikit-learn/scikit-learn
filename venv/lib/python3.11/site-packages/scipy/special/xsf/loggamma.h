/* Translated from Cython into C++ by SciPy developers in 2024.
 * Original header comment appears below.
 */

/* An implementation of the principal branch of the logarithm of
 * Gamma. Also contains implementations of Gamma and 1/Gamma which are
 * easily computed from log-Gamma.
 *
 * Author: Josh Wilson
 *
 * Distributed under the same license as Scipy.
 *
 * References
 * ----------
 * [1] Hare, "Computing the Principal Branch of log-Gamma",
 *     Journal of Algorithms, 1997.
 *
 * [2] Julia,
 *     https://github.com/JuliaLang/julia/blob/master/base/special/gamma.jl
 */

#pragma once

#include "cephes/gamma.h"
#include "cephes/rgamma.h"
#include "config.h"
#include "error.h"
#include "evalpoly.h"
#include "trig.h"
#include "zlog1.h"

namespace xsf {

namespace detail {
    constexpr double loggamma_SMALLX = 7;
    constexpr double loggamma_SMALLY = 7;
    constexpr double loggamma_HLOG2PI = 0.918938533204672742;      // log(2*pi)/2
    constexpr double loggamma_LOGPI = 1.1447298858494001741434262; // log(pi)
    constexpr double loggamma_TAYLOR_RADIUS = 0.2;

    XSF_HOST_DEVICE std::complex<double> loggamma_stirling(std::complex<double> z) {
        /* Stirling series for log-Gamma
         *
         * The coefficients are B[2*n]/(2*n*(2*n - 1)) where B[2*n] is the
         * (2*n)th Bernoulli number. See (1.1) in [1].
         */
        double coeffs[] = {-2.955065359477124183E-2,  6.4102564102564102564E-3, -1.9175269175269175269E-3,
                           8.4175084175084175084E-4,  -5.952380952380952381E-4, 7.9365079365079365079E-4,
                           -2.7777777777777777778E-3, 8.3333333333333333333E-2};
        std::complex<double> rz = 1.0 / z;
        std::complex<double> rzz = rz / z;

        return (z - 0.5) * std::log(z) - z + loggamma_HLOG2PI + rz * cevalpoly(coeffs, 7, rzz);
    }

    XSF_HOST_DEVICE std::complex<double> loggamma_recurrence(std::complex<double> z) {
        /* Backward recurrence relation.
         *
         * See Proposition 2.2 in [1] and the Julia implementation [2].
         *
         */
        int signflips = 0;
        int sb = 0;
        std::complex<double> shiftprod = z;

        z += 1.0;
        int nsb;
        while (z.real() <= loggamma_SMALLX) {
            shiftprod *= z;
            nsb = std::signbit(shiftprod.imag());
            signflips += nsb != 0 && sb == 0 ? 1 : 0;
            sb = nsb;
            z += 1.0;
        }
        return loggamma_stirling(z) - std::log(shiftprod) - signflips * 2 * M_PI * std::complex<double>(0, 1);
    }

    XSF_HOST_DEVICE std::complex<double> loggamma_taylor(std::complex<double> z) {
        /* Taylor series for log-Gamma around z = 1.
         *
         * It is
         *
         * loggamma(z + 1) = -gamma*z + zeta(2)*z**2/2 - zeta(3)*z**3/3 ...
         *
         * where gamma is the Euler-Mascheroni constant.
         */

        double coeffs[] = {
            -4.3478266053040259361E-2, 4.5454556293204669442E-2, -4.7619070330142227991E-2, 5.000004769810169364E-2,
            -5.2631679379616660734E-2, 5.5555767627403611102E-2, -5.8823978658684582339E-2, 6.2500955141213040742E-2,
            -6.6668705882420468033E-2, 7.1432946295361336059E-2, -7.6932516411352191473E-2, 8.3353840546109004025E-2,
            -9.0954017145829042233E-2, 1.0009945751278180853E-1, -1.1133426586956469049E-1, 1.2550966952474304242E-1,
            -1.4404989676884611812E-1, 1.6955717699740818995E-1, -2.0738555102867398527E-1, 2.7058080842778454788E-1,
            -4.0068563438653142847E-1, 8.2246703342411321824E-1, -5.7721566490153286061E-1};

        z -= 1.0;
        return z * cevalpoly(coeffs, 22, z);
    }
} // namespace detail

XSF_HOST_DEVICE inline double loggamma(double x) {
    if (x < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return cephes::lgam(x);
}

XSF_HOST_DEVICE inline float loggamma(float x) { return loggamma(static_cast<double>(x)); }

XSF_HOST_DEVICE inline std::complex<double> loggamma(std::complex<double> z) {
    // Compute the principal branch of log-Gamma

    if (std::isnan(z.real()) || std::isnan(z.imag())) {
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
    if (z.real() <= 0 and z == std::floor(z.real())) {
        set_error("loggamma", SF_ERROR_SINGULAR, NULL);
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
    if (z.real() > detail::loggamma_SMALLX || std::abs(z.imag()) > detail::loggamma_SMALLY) {
        return detail::loggamma_stirling(z);
    }
    if (std::abs(z - 1.0) < detail::loggamma_TAYLOR_RADIUS) {
        return detail::loggamma_taylor(z);
    }
    if (std::abs(z - 2.0) < detail::loggamma_TAYLOR_RADIUS) {
        // Recurrence relation and the Taylor series around 1.
        return detail::zlog1(z - 1.0) + detail::loggamma_taylor(z - 1.0);
    }
    if (z.real() < 0.1) {
        // Reflection formula; see Proposition 3.1 in [1]
        double tmp = std::copysign(2 * M_PI, z.imag()) * std::floor(0.5 * z.real() + 0.25);
        return std::complex<double>(detail::loggamma_LOGPI, tmp) - std::log(sinpi(z)) - loggamma(1.0 - z);
    }
    if (std::signbit(z.imag()) == 0) {
        // z.imag() >= 0 but is not -0.0
        return detail::loggamma_recurrence(z);
    }
    return std::conj(detail::loggamma_recurrence(std::conj(z)));
}

XSF_HOST_DEVICE inline std::complex<float> loggamma(std::complex<float> z) {
    return static_cast<std::complex<float>>(loggamma(static_cast<std::complex<double>>(z)));
}

XSF_HOST_DEVICE inline double rgamma(double z) { return cephes::rgamma(z); }

XSF_HOST_DEVICE inline float rgamma(float z) { return rgamma(static_cast<double>(z)); }

XSF_HOST_DEVICE inline std::complex<double> rgamma(std::complex<double> z) {
    // Compute 1/Gamma(z) using loggamma.
    if (z.real() <= 0 && z == std::floor(z.real())) {
        // Zeros at 0, -1, -2, ...
        return 0.0;
    }
    return std::exp(-loggamma(z));
}

XSF_HOST_DEVICE inline std::complex<float> rgamma(std::complex<float> z) {
    return static_cast<std::complex<float>>(rgamma(static_cast<std::complex<double>>(z)));
}

} // namespace xsf
