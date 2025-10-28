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

#include "cephes/sindg.h"
#include "cephes/tandg.h"
#include "cephes/trig.h"
#include "cephes/unity.h"
#include "config.h"
#include "evalpoly.h"

namespace xsf {

template <typename T>
XSF_HOST_DEVICE T sinpi(T x) {
    return cephes::sinpi(x);
}

template <typename T>
XSF_HOST_DEVICE std::complex<T> sinpi(std::complex<T> z) {
    T x = z.real();
    T piy = M_PI * z.imag();
    T abspiy = std::abs(piy);
    T sinpix = cephes::sinpi(x);
    T cospix = cephes::cospi(x);

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
    T exphpiy = std::exp(abspiy / 2);
    T coshfac;
    T sinhfac;
    if (exphpiy == std::numeric_limits<T>::infinity()) {
        if (sinpix == 0.0) {
            // Preserve the sign of zero.
            coshfac = std::copysign(0.0, sinpix);
        } else {
            coshfac = std::copysign(std::numeric_limits<T>::infinity(), sinpix);
        }
        if (cospix == 0.0) {
            // Preserve the sign of zero.
            sinhfac = std::copysign(0.0, cospix);
        } else {
            sinhfac = std::copysign(std::numeric_limits<T>::infinity(), cospix);
        }
        return {coshfac, sinhfac};
    }

    coshfac = 0.5 * sinpix * exphpiy;
    sinhfac = 0.5 * cospix * exphpiy;
    return {coshfac * exphpiy, sinhfac * exphpiy};
}

template <typename T>
XSF_HOST_DEVICE T cospi(T x) {
    return cephes::cospi(x);
}

template <typename T>
XSF_HOST_DEVICE std::complex<T> cospi(std::complex<T> z) {
    T x = z.real();
    T piy = M_PI * z.imag();
    T abspiy = std::abs(piy);
    T sinpix = cephes::sinpi(x);
    T cospix = cephes::cospi(x);

    if (abspiy < 700) {
        return {cospix * std::cosh(piy), -sinpix * std::sinh(piy)};
    }

    // See csinpi(z) for an idea of what's going on here.
    T exphpiy = std::exp(abspiy / 2);
    T coshfac;
    T sinhfac;
    if (exphpiy == std::numeric_limits<T>::infinity()) {
        if (sinpix == 0.0) {
            // Preserve the sign of zero.
            coshfac = std::copysign(0.0, cospix);
        } else {
            coshfac = std::copysign(std::numeric_limits<T>::infinity(), cospix);
        }
        if (cospix == 0.0) {
            // Preserve the sign of zero.
            sinhfac = std::copysign(0.0, sinpix);
        } else {
            sinhfac = std::copysign(std::numeric_limits<T>::infinity(), sinpix);
        }
        return {coshfac, sinhfac};
    }

    coshfac = 0.5 * cospix * exphpiy;
    sinhfac = 0.5 * sinpix * exphpiy;
    return {coshfac * exphpiy, sinhfac * exphpiy};
}

template <typename T>
XSF_HOST_DEVICE T sindg(T x) {
    return cephes::sindg(x);
}

template <>
XSF_HOST_DEVICE inline float sindg(float x) {
    return sindg(static_cast<double>(x));
}

template <typename T>
XSF_HOST_DEVICE T cosdg(T x) {
    return cephes::cosdg(x);
}

template <>
XSF_HOST_DEVICE inline float cosdg(float x) {
    return cosdg(static_cast<double>(x));
}

template <typename T>
XSF_HOST_DEVICE T tandg(T x) {
    return cephes::tandg(x);
}

template <>
XSF_HOST_DEVICE inline float tandg(float x) {
    return tandg(static_cast<double>(x));
}

template <typename T>
XSF_HOST_DEVICE T cotdg(T x) {
    return cephes::cotdg(x);
}

template <>
XSF_HOST_DEVICE inline float cotdg(float x) {
    return cotdg(static_cast<double>(x));
}

inline double radian(double d, double m, double s) { return cephes::radian(d, m, s); }

inline float radian(float d, float m, float s) {
    return radian(static_cast<double>(d), static_cast<double>(m), static_cast<double>(s));
}

inline double cosm1(double x) { return cephes::cosm1(x); }

inline float cosm1(float x) { return cosm1(static_cast<double>(x)); }

} // namespace xsf
