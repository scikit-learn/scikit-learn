/* Translated from Cython into C++ by SciPy developers in 2024.
 * Original header comment appears below.
 */

/* An implementation of the digamma function for complex arguments.
 *
 * Author: Josh Wilson
 *
 * Distributed under the same license as Scipy.
 *
 * Sources:
 * [1] "The Digital Library of Mathematical Functions", dlmf.nist.gov
 *
 * [2] mpmath (version 0.19), http://mpmath.org
 */

#pragma once

#include "cephes/psi.h"
#include "cephes/zeta.h"
#include "config.h"
#include "error.h"
#include "trig.h"

namespace special {
namespace detail {
    // All of the following were computed with mpmath
    // Location of the positive root
    constexpr double digamma_posroot = 1.4616321449683623;
    // Value of the positive root
    constexpr double digamma_posrootval = -9.2412655217294275e-17;
    // Location of the negative root
    constexpr double digamma_negroot = -0.504083008264455409;
    // Value of the negative root
    constexpr double digamma_negrootval = 7.2897639029768949e-17;

    template <typename T>
    SPECFUN_HOST_DEVICE T digamma_zeta_series(T z, double root, double rootval) {
        T res = rootval;
        T coeff = -1.0;

        z = z - root;
        T term;
        for (int n = 1; n < 100; n++) {
            coeff *= -z;
            term = coeff * cephes::zeta(n + 1, root);
            res += term;
            if (std::abs(term) < std::numeric_limits<double>::epsilon() * std::abs(res)) {
                break;
            }
        }
        return res;
    }

    SPECFUN_HOST_DEVICE inline std::complex<double> digamma_forward_recurrence(std::complex<double> z,
                                                                               std::complex<double> psiz, int n) {
        /* Compute digamma(z + n) using digamma(z) using the recurrence
         * relation
         *
         * digamma(z + 1) = digamma(z) + 1/z.
         *
         * See https://dlmf.nist.gov/5.5#E2 */
        std::complex<double> res = psiz;

        for (int k = 0; k < n; k++) {
            res += 1.0 / (z + static_cast<double>(k));
        }
        return res;
    }

    SPECFUN_HOST_DEVICE inline std::complex<double> digamma_backward_recurrence(std::complex<double> z,
                                                                                std::complex<double> psiz, int n) {
        /* Compute digamma(z - n) using digamma(z) and a recurrence relation. */
        std::complex<double> res = psiz;

        for (int k = 1; k < n + 1; k++) {
            res -= 1.0 / (z - static_cast<double>(k));
        }
        return res;
    }

    SPECFUN_HOST_DEVICE inline std::complex<double> digamma_asymptotic_series(std::complex<double> z) {
        /* Evaluate digamma using an asymptotic series. See
         *
         * https://dlmf.nist.gov/5.11#E2 */
        double bernoulli2k[] = {
            0.166666666666666667,  -0.0333333333333333333, 0.0238095238095238095, -0.0333333333333333333,
            0.0757575757575757576, -0.253113553113553114,  1.16666666666666667,   -7.09215686274509804,
            54.9711779448621554,   -529.124242424242424,   6192.12318840579710,   -86580.2531135531136,
            1425517.16666666667,   -27298231.0678160920,   601580873.900642368,   -15116315767.0921569};
        std::complex<double> rzz = 1.0 / z / z;
        std::complex<double> zfac = 1.0;
        std::complex<double> term;
        std::complex<double> res;

        if (!(std::isfinite(z.real()) && std::isfinite(z.imag()))) {
            /* Check for infinity (or nan) and return early.
             * Result of division by complex infinity is implementation dependent.
             * and has been observed to vary between C++ stdlib and CUDA stdlib.
             */
            return std::log(z);
        }

        res = std::log(z) - 0.5 / z;

        for (int k = 1; k < 17; k++) {
            zfac *= rzz;
            term = -bernoulli2k[k - 1] * zfac / (2 * static_cast<double>(k));
            res += term;
            if (std::abs(term) < std::numeric_limits<double>::epsilon() * std::abs(res)) {
                break;
            }
        }
        return res;
    }

} // namespace detail

SPECFUN_HOST_DEVICE inline double digamma(double z) {
    /* Wrap Cephes' psi to take advantage of the series expansion around
     * the smallest negative zero.
     */
    if (std::abs(z - detail::digamma_negroot) < 0.3) {
        return detail::digamma_zeta_series(z, detail::digamma_negroot, detail::digamma_negrootval);
    }
    return cephes::psi(z);
}

SPECFUN_HOST_DEVICE inline float digamma(float z) { return static_cast<float>(digamma(static_cast<double>(z))); }

SPECFUN_HOST_DEVICE inline std::complex<double> digamma(std::complex<double> z) {
    /*
     * Compute the digamma function for complex arguments. The strategy
     * is:
     *
     * - Around the two zeros closest to the origin (posroot and negroot)
     * use a Taylor series with precomputed zero order coefficient.
     * - If close to the origin, use a recurrence relation to step away
     * from the origin.
     * - If close to the negative real axis, use the reflection formula
     * to move to the right halfplane.
     * - If |z| is large (> 16), use the asymptotic series.
     * - If |z| is small, use a recurrence relation to make |z| large
     * enough to use the asymptotic series.
     */
    double absz = std::abs(z);
    std::complex<double> res = 0;
    /* Use the asymptotic series for z away from the negative real axis
     * with abs(z) > smallabsz. */
    int smallabsz = 16;
    /* Use the reflection principle for z with z.real < 0 that are within
     * smallimag of the negative real axis.
     * int smallimag = 6  # unused below except in a comment */

    if (z.real() <= 0.0 && std::ceil(z.real()) == z) {
        // Poles
        set_error("digamma", SF_ERROR_SINGULAR, NULL);
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
    if (std::abs(z - detail::digamma_negroot) < 0.3) {
        // First negative root.
        return detail::digamma_zeta_series(z, detail::digamma_negroot, detail::digamma_negrootval);
    }

    if (z.real() < 0 and std::abs(z.imag()) < smallabsz) {
        /* Reflection formula for digamma. See
         *
         *https://dlmf.nist.gov/5.5#E4
         */
        res = -M_PI * cospi(z) / sinpi(z);
        z = 1.0 - z;
        absz = std::abs(z);
    }

    if (absz < 0.5) {
        /* Use one step of the recurrence relation to step away from
         * the pole. */
        res = -1.0 / z;
        z += 1.0;
        absz = std::abs(z);
    }

    if (std::abs(z - detail::digamma_posroot) < 0.5) {
        res += detail::digamma_zeta_series(z, detail::digamma_posroot, detail::digamma_posrootval);
    } else if (absz > smallabsz) {
        res += detail::digamma_asymptotic_series(z);
    } else if (z.real() >= 0.0) {
        double n = std::trunc(smallabsz - absz) + 1;
        std::complex<double> init = detail::digamma_asymptotic_series(z + n);
        res += detail::digamma_backward_recurrence(z + n, init, n);
    } else {
        // z.real() < 0, absz < smallabsz, and z.imag() > smallimag
        double n = std::trunc(smallabsz - absz) - 1;
        std::complex<double> init = detail::digamma_asymptotic_series(z - n);
        res += detail::digamma_forward_recurrence(z - n, init, n);
    }
    return res;
}

SPECFUN_HOST_DEVICE inline std::complex<float> digamma(std::complex<float> z) {
    return static_cast<std::complex<float>>(digamma(static_cast<std::complex<double>>(z)));
}

} // namespace special
