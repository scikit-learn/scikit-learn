/* The functions exp1, expi below are based on translations of the Fortran code
 * by Shanjie Zhang and Jianming Jin from the book
 *
 *  Shanjie Zhang, Jianming Jin,
 *  Computation of Special Functions,
 *  Wiley, 1996,
 *  ISBN: 0-471-11963-6,
 *  LC: QA351.C45.
 */

#pragma once

#include "config.h"
#include "error.h"

#include "cephes/const.h"


namespace xsf {


XSF_HOST_DEVICE inline double exp1(double x) {
    // ============================================
    // Purpose: Compute exponential integral E1(x)
    // Input :  x  --- Argument of E1(x)
    // Output:  E1 --- E1(x)  ( x > 0 )
    // ============================================
    int k, m;
    double e1, r, t, t0;
    constexpr double ga = cephes::detail::SCIPY_EULER;

    if (x == 0.0) {
	return std::numeric_limits<double>::infinity();
    }
    if (x <= 1.0) {
        e1 = 1.0;
        r = 1.0;
        for (k = 1; k < 26; k++) {
            r = -r*k*x/std::pow(k+1.0, 2);
            e1 += r;
            if (std::abs(r) <= std::abs(e1)*1e-15) { break; }
        }
        return -ga - std::log(x) + x*e1;
    }
    m = 20 + (int)(80.0/x);
    t0 = 0.0;
    for (k = m; k > 0; k--) {
	t0 = k / (1.0 + k / (x+t0));
    }
    t = 1.0 / (x + t0);
    return std::exp(-x)*t;
}

XSF_HOST_DEVICE inline float exp1(float x) { return exp1(static_cast<double>(x)); }

XSF_HOST_DEVICE inline std::complex<double> exp1(std::complex<double> z) {
    // ====================================================
    // Purpose: Compute complex exponential integral E1(z)
    // Input :  z   --- Argument of E1(z)
    // Output:  CE1 --- E1(z)
    // ====================================================
    constexpr double el = cephes::detail::SCIPY_EULER;
    int k;
    std::complex<double> ce1, cr, zc, zd, zdc;
    double x = z.real();
    double a0 = std::abs(z);
    // Continued fraction converges slowly near negative real axis,
    // so use power series in a wedge around it until radius 40.0
    double xt = -2.0*std::abs(z.imag());

    if (a0 == 0.0) { return std::numeric_limits<double>::infinity(); }
    if ((a0 < 5.0) || ((x < xt) && (a0 < 40.0))) {
        // Power series
        ce1 = 1.0;
        cr = 1.0;
        for (k = 1; k < 501; k++) {
            cr = -cr*z*static_cast<double>(k / std::pow(k + 1, 2));
            ce1 += cr;
            if (std::abs(cr) < std::abs(ce1)*1e-15) { break; }
        }
        if ((x <= 0.0) && (z.imag() == 0.0)) {
            //Careful on the branch cut -- use the sign of the imaginary part
            // to get the right sign on the factor if pi.
            ce1 = -el - std::log(-z) + z*ce1 - std::copysign(M_PI, z.imag())*std::complex<double>(0.0, 1.0);
        } else {
            ce1 = -el - std::log(z) + z*ce1;
        }
    } else {
        // Continued fraction https://dlmf.nist.gov/6.9
        //                  1     1     1     2     2     3     3
        // E1 = exp(-z) * ----- ----- ----- ----- ----- ----- ----- ...
        //                Z +   1 +   Z +   1 +   Z +   1 +   Z +
        zc = 0.0;
        zd = static_cast<double>(1) / z;
        zdc = zd;
        zc += zdc;
        for (k = 1; k < 501; k++) {
            zd = static_cast<double>(1) / (zd*static_cast<double>(k) + static_cast<double>(1));
            zdc *= (zd - static_cast<double>(1));
            zc += zdc;

            zd = static_cast<double>(1) / (zd*static_cast<double>(k) + z);
            zdc *= (z*zd - static_cast<double>(1));
            zc += zdc;
            if ((std::abs(zdc) <= std::abs(zc)*1e-15) && (k > 20)) { break; }
        }
        ce1 = std::exp(-z)*zc;
        if ((x <= 0.0) && (z.imag() == 0.0)) {
            ce1 -= M_PI*std::complex<double>(0.0, 1.0);
        }
    }
    return ce1;
}

XSF_HOST_DEVICE inline std::complex<float> exp1(std::complex<float> z) {
    return static_cast<std::complex<float>>(exp1(static_cast<std::complex<double>>(z)));
}

XSF_HOST_DEVICE inline double expi(double x) {
    // ============================================
    // Purpose: Compute exponential integral Ei(x)
    // Input :  x  --- Argument of Ei(x)
    // Output:  EI --- Ei(x)
    // ============================================

    constexpr double ga = cephes::detail::SCIPY_EULER;
    double ei, r;

    if (x == 0.0) {
        ei = -std::numeric_limits<double>::infinity();
    } else if (x < 0) {
        ei = -exp1(-x);
    } else if (std::abs(x) <= 40.0) {
        // Power series around x=0
        ei = 1.0;
        r = 1.0;

        for (int k = 1; k <= 100; k++) {
            r = r * k * x / ((k + 1.0) * (k + 1.0));
            ei += r;
            if (std::abs(r / ei) <= 1.0e-15) { break; }
        }
        ei = ga + std::log(x) + x * ei;
    } else {
        // Asymptotic expansion (the series is not convergent)
        ei = 1.0;
        r = 1.0;
        for (int k = 1; k <= 20; k++) {
            r = r * k / x;
            ei += r;
        }
        ei = std::exp(x) / x * ei;
    }
    return ei;
}

XSF_HOST_DEVICE inline float expi(float x) { return expi(static_cast<double>(x)); }
    
std::complex<double> expi(std::complex<double> z) {
    // ============================================
    // Purpose: Compute exponential integral Ei(x)
    // Input :  x  --- Complex argument of Ei(x)
    // Output:  EI --- Ei(x)
    // ============================================

    std::complex<double> cei;
    cei = - exp1(-z);
    if (z.imag() > 0.0) {
        cei += std::complex<double>(0.0, M_PI);
    } else if (z.imag() < 0.0 ) {
        cei -= std::complex<double>(0.0, M_PI);
    } else {
        if (z.real() > 0.0) {
            cei += std::complex<double>(0.0, copysign(M_PI, z.imag()));
        }
    }
    return cei;
}


XSF_HOST_DEVICE inline std::complex<float> expi(std::complex<float> z) {
    return static_cast<std::complex<float>>(expi(static_cast<std::complex<double>>(z)));
}

namespace detail {

    //
    // Compute a factor of the exponential integral E1.
    // This is used in scaled_exp1(x) for moderate values of x.
    //
    // The function uses the continued fraction expansion given in equation 5.1.22
    // of Abramowitz & Stegun, "Handbook of Mathematical Functions".
    // For n=1, this is
    //
    //    E1(x) = exp(-x)*C(x)
    //
    // where C(x), expressed in the notation used in A&S, is the continued fraction
    //
    //            1    1    1    2    2    3    3
    //    C(x) = ---  ---  ---  ---  ---  ---  ---  ...
    //           x +  1 +  x +  1 +  x +  1 +  x +
    //
    // Here, we pull a factor of 1/z out of C(x), so
    //
    //    E1(x) = (exp(-x)/x)*F(x)
    //
    // and a bit of algebra gives the continued fraction expansion of F(x) to be
    //
    //            1    1    1    2    2    3    3
    //    F(x) = ---  ---  ---  ---  ---  ---  ---  ...
    //           1 +  x +  1 +  x +  1 +  x +  1 +
    //
    XSF_HOST_DEVICE inline double expint1_factor_cont_frac(double x) {
        // The number of terms to use in the truncated continued fraction
        // depends on x.  Larger values of x require fewer terms.
        int m = 20 + (int) (80.0 / x);
        double t0 = 0.0;
        for (int k = m; k > 0; --k) {
            t0 = k / (x + k / (1 + t0));
        }
        return 1 / (1 + t0);
    }

} // namespace detail

//
// Scaled version  of the exponential integral E_1(x).
//
// Factor E_1(x) as
//
//    E_1(x) = exp(-x)/x * F(x)
//
// This function computes F(x).
//
// F(x) has the properties:
//  * F(0) = 0
//  * F is increasing on [0, inf)
//  * lim_{x->inf} F(x) = 1.
//
XSF_HOST_DEVICE inline double scaled_exp1(double x) {
    if (x < 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (x == 0) {
        return 0.0;
    }

    if (x <= 1) {
        // For small x, the naive implementation is sufficiently accurate.
        return x * std::exp(x) * exp1(x);
    }

    if (x <= 1250) {
        // For moderate x, use the continued fraction expansion.
        return detail::expint1_factor_cont_frac(x);
    }

    // For large x, use the asymptotic expansion.  This is equation 5.1.51
    // from Abramowitz & Stegun, "Handbook of Mathematical Functions".
    return 1 + (-1 + (2 + (-6 + (24 - 120 / x) / x) / x) / x) / x;
}

XSF_HOST_DEVICE inline float scaled_exp1(float x) { return scaled_exp1(static_cast<double>(x)); }

} // namespace xsf
