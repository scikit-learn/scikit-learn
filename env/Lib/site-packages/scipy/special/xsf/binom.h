/* Translated from Cython into C++ by SciPy developers in 2024.
 *
 * Original authors: Pauli Virtanen, Eric Moore
 */

// Binomial coefficient

#pragma once

#include "config.h"

#include "cephes/beta.h"
#include "cephes/gamma.h"

namespace xsf {

XSF_HOST_DEVICE inline double binom(double n, double k) {
    double kx, nx, num, den, dk, sgn;

    if (n < 0) {
        nx = std::floor(n);
        if (n == nx) {
            // Undefined
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    kx = std::floor(k);
    if (k == kx && (std::abs(n) > 1E-8 || n == 0)) {
        /* Integer case: use multiplication formula for less rounding
         * error for cases where the result is an integer.
         *
         * This cannot be used for small nonzero n due to loss of
         * precision. */
        nx = std::floor(n);
        if (nx == n && kx > nx / 2 && nx > 0) {
            // Reduce kx by symmetry
            kx = nx - kx;
        }

        if (kx >= 0 && kx < 20) {
            num = 1.0;
            den = 1.0;
            for (int i = 1; i < 1 + static_cast<int>(kx); i++) {
                num *= i + n - kx;
                den *= i;
                if (std::abs(num) > 1E50) {
                    num /= den;
                    den = 1.0;
                }
            }
            return num / den;
        }
    }

    // general case
    if (n >= 1E10 * k and k > 0) {
        // avoid under/overflows intermediate results
        return std::exp(-cephes::lbeta(1 + n - k, 1 + k) - std::log(n + 1));
    }
    if (k > 1E8 * std::abs(n)) {
        // avoid loss of precision
        num = cephes::Gamma(1 + n) / std::abs(k) + cephes::Gamma(1 + n) * n / (2 * k * k); // + ...
        num /= M_PI * std::pow(std::abs(k), n);
        if (k > 0) {
            kx = std::floor(k);
            if (static_cast<int>(kx) == kx) {
                dk = k - kx;
                sgn = (static_cast<int>(kx) % 2 == 0) ? 1 : -1;
            } else {
                dk = k;
                sgn = 1;
            }
            return num * std::sin((dk - n) * M_PI) * sgn;
        }
        kx = std::floor(k);
        if (static_cast<int>(kx) == kx) {
            return 0;
        }
        return num * std::sin(k * M_PI);
    }
    return 1 / (n + 1) / cephes::beta(1 + n - k, 1 + k);
}

XSF_HOST_DEVICE inline float binom(float n, float k) {
    return binom(static_cast<double>(n), static_cast<double>(k));
}

} // namespace xsf
