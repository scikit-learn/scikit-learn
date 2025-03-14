/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     beta.c
 *
 *     Beta function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, y, beta();
 *
 * y = beta( a, b );
 *
 *
 *
 * DESCRIPTION:
 *
 *                   -     -
 *                  | (a) | (b)
 * beta( a, b )  =  -----------.
 *                     -
 *                    | (a+b)
 *
 * For large arguments the logarithm of the function is
 * evaluated using lgam(), then exponentiated.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,30       30000       8.1e-14     1.1e-14
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * beta overflow    log(beta) > MAXLOG       0.0
 *                  a or b <0 integer        0.0
 *
 */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
#pragma once

#include "../config.h"
#include "const.h"
#include "gamma.h"
#include "rgamma.h"

namespace xsf {
namespace cephes {

    XSF_HOST_DEVICE double beta(double, double);
    XSF_HOST_DEVICE double lbeta(double, double);

    namespace detail {
        constexpr double beta_ASYMP_FACTOR = 1e6;

        /*
         * Asymptotic expansion for  ln(|B(a, b)|) for a > ASYMP_FACTOR*max(|b|, 1).
         */
        XSF_HOST_DEVICE inline double lbeta_asymp(double a, double b, int *sgn) {
            double r = lgam_sgn(b, sgn);
            r -= b * std::log(a);

            r += b * (1 - b) / (2 * a);
            r += b * (1 - b) * (1 - 2 * b) / (12 * a * a);
            r += -b * b * (1 - b) * (1 - b) / (12 * a * a * a);

            return r;
        }

        /*
         * Special case for a negative integer argument
         */

        XSF_HOST_DEVICE inline double beta_negint(int a, double b) {
            int sgn;
            if (b == static_cast<int>(b) && 1 - a - b > 0) {
                sgn = (static_cast<int>(b) % 2 == 0) ? 1 : -1;
                return sgn * xsf::cephes::beta(1 - a - b, b);
            } else {
                set_error("lbeta", SF_ERROR_OVERFLOW, NULL);
                return std::numeric_limits<double>::infinity();
            }
        }

        XSF_HOST_DEVICE inline double lbeta_negint(int a, double b) {
            double r;
            if (b == static_cast<int>(b) && 1 - a - b > 0) {
                r = xsf::cephes::lbeta(1 - a - b, b);
                return r;
            } else {
                set_error("lbeta", SF_ERROR_OVERFLOW, NULL);
                return std::numeric_limits<double>::infinity();
            }
        }
    } // namespace detail

    XSF_HOST_DEVICE inline double beta(double a, double b) {
        double y;
        int sign = 1;

        if (a <= 0.0) {
            if (a == std::floor(a)) {
                if (a == static_cast<int>(a)) {
                    return detail::beta_negint(static_cast<int>(a), b);
                } else {
                    goto overflow;
                }
            }
        }

        if (b <= 0.0) {
            if (b == std::floor(b)) {
                if (b == static_cast<int>(b)) {
                    return detail::beta_negint(static_cast<int>(b), a);
                } else {
                    goto overflow;
                }
            }
        }

        if (std::abs(a) < std::abs(b)) {
            y = a;
            a = b;
            b = y;
        }

        if (std::abs(a) > detail::beta_ASYMP_FACTOR * std::abs(b) && a > detail::beta_ASYMP_FACTOR) {
            /* Avoid loss of precision in lgam(a + b) - lgam(a) */
            y = detail::lbeta_asymp(a, b, &sign);
            return sign * std::exp(y);
        }

        y = a + b;
        if (std::abs(y) > detail::MAXGAM || std::abs(a) > detail::MAXGAM || std::abs(b) > detail::MAXGAM) {
            int sgngam;
            y = detail::lgam_sgn(y, &sgngam);
            sign *= sgngam; /* keep track of the sign */
            y = detail::lgam_sgn(b, &sgngam) - y;
            sign *= sgngam;
            y = detail::lgam_sgn(a, &sgngam) + y;
            sign *= sgngam;
            if (y > detail::MAXLOG) {
                goto overflow;
            }
            return (sign * std::exp(y));
        }

        y = rgamma(y);
        a = Gamma(a);
        b = Gamma(b);
        if (std::isinf(y)) {
            goto overflow;
	}

        if (std::abs(std::abs(a*y) - 1.0) > std::abs(std::abs(b*y) - 1.0)) {
            y = b * y;
            y *= a;
        } else {
            y = a * y;
            y *= b;
        }

        return (y);

    overflow:
        set_error("beta", SF_ERROR_OVERFLOW, NULL);
        return (sign * std::numeric_limits<double>::infinity());
    }

    /* Natural log of |beta|. */

    XSF_HOST_DEVICE inline double lbeta(double a, double b) {
        double y;
        int sign;

        sign = 1;

        if (a <= 0.0) {
            if (a == std::floor(a)) {
                if (a == static_cast<int>(a)) {
                    return detail::lbeta_negint(static_cast<int>(a), b);
                } else {
                    goto over;
                }
            }
        }

        if (b <= 0.0) {
            if (b == std::floor(b)) {
                if (b == static_cast<int>(b)) {
                    return detail::lbeta_negint(static_cast<int>(b), a);
                } else {
                    goto over;
                }
            }
        }

        if (std::abs(a) < std::abs(b)) {
            y = a;
            a = b;
            b = y;
        }

        if (std::abs(a) > detail::beta_ASYMP_FACTOR * std::abs(b) && a > detail::beta_ASYMP_FACTOR) {
            /* Avoid loss of precision in lgam(a + b) - lgam(a) */
            y = detail::lbeta_asymp(a, b, &sign);
            return y;
        }

        y = a + b;
        if (std::abs(y) > detail::MAXGAM || std::abs(a) > detail::MAXGAM || std::abs(b) > detail::MAXGAM) {
            int sgngam;
            y = detail::lgam_sgn(y, &sgngam);
            sign *= sgngam; /* keep track of the sign */
            y = detail::lgam_sgn(b, &sgngam) - y;
            sign *= sgngam;
            y = detail::lgam_sgn(a, &sgngam) + y;
            sign *= sgngam;
            return (y);
        }

        y = rgamma(y);
        a = Gamma(a);
        b = Gamma(b);
        if (std::isinf(y)) {
        over:
            set_error("lbeta", SF_ERROR_OVERFLOW, NULL);
            return (sign * std::numeric_limits<double>::infinity());
        }

        if (std::abs(std::abs(a*y) - 1.0) > std::abs(std::abs(b*y) - 1.0)) {
            y = b * y;
            y *= a;
        } else {
            y = a * y;
            y *= b;
        }

        if (y < 0) {
            y = -y;
        }

        return (std::log(y));
    }
} // namespace cephes
} // namespace xsf
