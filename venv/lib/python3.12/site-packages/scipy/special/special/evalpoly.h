/* Translated from Cython into C++ by SciPy developers in 2024.
 *
 * Original author: Josh Wilson, 2016.
 */

/* Evaluate polynomials.
 *
 * All of the coefficients are stored in reverse order, i.e. if the
 * polynomial is
 *
 *     u_n x^n + u_{n - 1} x^{n - 1} + ... + u_0,
 *
 * then coeffs[0] = u_n, coeffs[1] = u_{n - 1}, ..., coeffs[n] = u_0.
 *
 * References
 * ----------
 * [1] Knuth, "The Art of Computer Programming, Volume II"
 */

#pragma once

#include "config.h"

namespace special {

SPECFUN_HOST_DEVICE inline std::complex<double> cevalpoly(const double *coeffs, int degree, std::complex<double> z) {
    /* Evaluate a polynomial with real coefficients at a complex point.
     *
     * Uses equation (3) in section 4.6.4 of [1]. Note that it is more
     * efficient than Horner's method.
     */
    double a = coeffs[0];
    double b = coeffs[1];
    double r = 2 * z.real();
    double s = std::norm(z);
    double tmp;

    for (int j = 2; j < degree + 1; j++) {
        tmp = b;
        b = std::fma(-s, a, coeffs[j]);
        a = std::fma(r, a, tmp);
    }

    return z * a + b;
}

} // namespace special
