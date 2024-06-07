/* Translated from Cython into C++ by SciPy developers in 2023.
 * Original header with Copyright information appears below.
 */

/* Implementation of the Lambert W function [1]. Based on MPMath
 *  Implementation [2], and documentation [3].
 *
 * Copyright: Yosef Meller, 2009
 * Author email: mellerf@netvision.net.il
 *
 * Distributed under the same license as SciPy
 *
 *
 * References:
 * [1] On the Lambert W function, Adv. Comp. Math. 5 (1996) 329-359,
 *     available online: https://web.archive.org/web/20230123211413/https://cs.uwaterloo.ca/research/tr/1993/03/W.pdf
 * [2] mpmath source code,
 https://github.com/mpmath/mpmath/blob/c5939823669e1bcce151d89261b802fe0d8978b4/mpmath/functions/functions.py#L435-L461
 * [3]
 https://web.archive.org/web/20230504171447/https://mpmath.org/doc/current/functions/powers.html#lambert-w-function
 *

 * TODO: use a series expansion when extremely close to the branch point
 * at `-1/e` and make sure that the proper branch is chosen there.
 */

#pragma once

#include "config.h"
#include "error.h"
#include "evalpoly.h"

namespace special {
constexpr double EXPN1 = 0.36787944117144232159553; // exp(-1)
constexpr double OMEGA = 0.56714329040978387299997; // W(1, 0)

namespace detail {
    SPECFUN_HOST_DEVICE inline std::complex<double> lambertw_branchpt(std::complex<double> z) {
        // Series for W(z, 0) around the branch point; see 4.22 in [1].
        double coeffs[] = {-1.0 / 3.0, 1.0, -1.0};
        std::complex<double> p = std::sqrt(2.0 * (M_E * z + 1.0));

        return cevalpoly(coeffs, 2, p);
    }

    SPECFUN_HOST_DEVICE inline std::complex<double> lambertw_pade0(std::complex<double> z) {
        // (3, 2) Pade approximation for W(z, 0) around 0.
        double num[] = {12.85106382978723404255, 12.34042553191489361902, 1.0};
        double denom[] = {32.53191489361702127660, 14.34042553191489361702, 1.0};

        /* This only gets evaluated close to 0, so we don't need a more
         * careful algorithm that avoids overflow in the numerator for
         * large z. */
        return z * cevalpoly(num, 2, z) / cevalpoly(denom, 2, z);
    }

    SPECFUN_HOST_DEVICE inline std::complex<double> lambertw_asy(std::complex<double> z, long k) {
        /* Compute the W function using the first two terms of the
         * asymptotic series. See 4.20 in [1].
         */
        std::complex<double> w = std::log(z) + 2.0 * M_PI * k * std::complex<double>(0, 1);
        return w - std::log(w);
    }

} // namespace detail

SPECFUN_HOST_DEVICE inline std::complex<double> lambertw(std::complex<double> z, long k, double tol) {
    double absz;
    std::complex<double> w;
    std::complex<double> ew, wew, wewz, wn;

    if (std::isnan(z.real()) || std::isnan(z.imag())) {
        return z;
    }
    if (z.real() == std::numeric_limits<double>::infinity()) {
        return z + 2.0 * M_PI * k * std::complex<double>(0, 1);
    }
    if (z.real() == -std::numeric_limits<double>::infinity()) {
        return -z + (2.0 * M_PI * k + M_PI) * std::complex<double>(0, 1);
    }
    if (z == 0.0) {
        if (k == 0) {
            return z;
        }
        set_error("lambertw", SF_ERROR_SINGULAR, NULL);
        return -std::numeric_limits<double>::infinity();
    }
    if (z == 1.0 && k == 0) {
        // Split out this case because the asymptotic series blows up
        return OMEGA;
    }

    absz = std::abs(z);
    // Get an initial guess for Halley's method
    if (k == 0) {
        if (std::abs(z + EXPN1) < 0.3) {
            w = detail::lambertw_branchpt(z);
        } else if (-1.0 < z.real() && z.real() < 1.5 && std::abs(z.imag()) < 1.0 &&
                   -2.5 * std::abs(z.imag()) - 0.2 < z.real()) {
            /* Empirically determined decision boundary where the Pade
             * approximation is more accurate. */
            w = detail::lambertw_pade0(z);
        } else {
            w = detail::lambertw_asy(z, k);
        }
    } else if (k == -1) {
        if (absz <= EXPN1 && z.imag() == 0.0 && z.real() < 0.0) {
            w = std::log(-z.real());
        } else {
            w = detail::lambertw_asy(z, k);
        }
    } else {
        w = detail::lambertw_asy(z, k);
    }

    // Halley's method; see 5.9 in [1]
    if (w.real() >= 0) {
        // Rearrange the formula to avoid overflow in exp
        for (int i = 0; i < 100; i++) {
            ew = std::exp(-w);
            wewz = w - z * ew;
            wn = w - wewz / (w + 1.0 - (w + 2.0) * wewz / (2.0 * w + 2.0));
            if (std::abs(wn - w) <= tol * std::abs(wn)) {
                return wn;
            }
            w = wn;
        }
    } else {
        for (int i = 0; i < 100; i++) {
            ew = std::exp(w);
            wew = w * ew;
            wewz = wew - z;
            wn = w - wewz / (wew + ew - (w + 2.0) * wewz / (2.0 * w + 2.0));
            if (std::abs(wn - w) <= tol * std::abs(wn)) {
                return wn;
            }
            w = wn;
        }
    }

    set_error("lambertw", SF_ERROR_SLOW, "iteration failed to converge: %g + %gj", z.real(), z.imag());
    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
}

} // namespace special
