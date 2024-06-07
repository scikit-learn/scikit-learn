/* Translated from Cython into C++ by SciPy developers in 2023.
 *
 * Original author: Josh Wilson, 2016.
 */

/* Implement sin(pi*z) and cos(pi*z) for complex z. Since the periods
 * of these functions are integral (and thus better representable in
 * floating point), it's possible to compute them with greater accuracy
 * than sin(z), cos(z).
 */

#pragma once

#include "cephes/trig.h"
#include "config.h"
#include "evalpoly.h"

namespace special {

SPECFUN_HOST_DEVICE inline std::complex<double> sinpi(std::complex<double> z) {
    double x = z.real();
    double piy = M_PI * z.imag();
    double abspiy = std::abs(piy);
    double sinpix = cephes::sinpi(x);
    double cospix = cephes::cospi(x);

    if (abspiy < 700) {
        return {sinpix * std::cosh(piy), cospix * std::sinh(piy)};
    }

    /* Have to be careful--sinh/cosh could overflow while cos/sin are small.
     * At this large of values
     *
     * cosh(y) ~ exp(y)/2
     * sinh(y) ~ sgn(y)*exp(y)/2
     *
     * so we can compute exp(y/2), scale by the right factor of sin/cos
     * and then multiply by exp(y/2) to avoid overflow. */
    double exphpiy = std::exp(abspiy / 2);
    double coshfac;
    double sinhfac;
    if (exphpiy == std::numeric_limits<double>::infinity()) {
        if (sinpix == 0.0) {
            // Preserve the sign of zero.
            coshfac = std::copysign(0.0, sinpix);
        } else {
            coshfac = std::copysign(std::numeric_limits<double>::infinity(), sinpix);
        }
        if (cospix == 0.0) {
            // Preserve the sign of zero.
            sinhfac = std::copysign(0.0, cospix);
        } else {
            sinhfac = std::copysign(std::numeric_limits<double>::infinity(), cospix);
        }
        return {coshfac, sinhfac};
    }

    coshfac = 0.5 * sinpix * exphpiy;
    sinhfac = 0.5 * cospix * exphpiy;
    return {coshfac * exphpiy, sinhfac * exphpiy};
}

SPECFUN_HOST_DEVICE inline std::complex<double> cospi(std::complex<double> z) {
    double x = z.real();
    double piy = M_PI * z.imag();
    double abspiy = std::abs(piy);
    double sinpix = cephes::sinpi(x);
    double cospix = cephes::cospi(x);

    if (abspiy < 700) {
        return {cospix * std::cosh(piy), -sinpix * std::sinh(piy)};
    }

    // See csinpi(z) for an idea of what's going on here.
    double exphpiy = std::exp(abspiy / 2);
    double coshfac;
    double sinhfac;
    if (exphpiy == std::numeric_limits<double>::infinity()) {
        if (sinpix == 0.0) {
            // Preserve the sign of zero.
            coshfac = std::copysign(0.0, cospix);
        } else {
            coshfac = std::copysign(std::numeric_limits<double>::infinity(), cospix);
        }
        if (cospix == 0.0) {
            // Preserve the sign of zero.
            sinhfac = std::copysign(0.0, sinpix);
        } else {
            sinhfac = std::copysign(std::numeric_limits<double>::infinity(), sinpix);
        }
        return {coshfac, sinhfac};
    }

    coshfac = 0.5 * cospix * exphpiy;
    sinhfac = 0.5 * sinpix * exphpiy;
    return {coshfac * exphpiy, sinhfac * exphpiy};
}

} // namespace special
