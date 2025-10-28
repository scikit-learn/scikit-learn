/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     ellpk.c
 *
 *     Complete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double m1, y, ellpk();
 *
 * y = ellpk( m1 );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *            pi/2
 *             -
 *            | |
 *            |           dt
 * K(m)  =    |    ------------------
 *            |                   2
 *          | |    sqrt( 1 - m sin t )
 *           -
 *            0
 *
 * where m = 1 - m1, using the approximation
 *
 *     P(x)  -  log x Q(x).
 *
 * The argument m1 is used internally rather than m so that the logarithmic
 * singularity at m = 1 will be shifted to the origin; this
 * preserves maximum accuracy.
 *
 * K(0) = pi/2.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,1        30000       2.5e-16     6.8e-17
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * ellpk domain       x<0, x>1           0.0
 *
 */

/*                                                     ellpk.c */

/*
 * Cephes Math Library, Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
#pragma once

#include "../config.h"
#include "../error.h"
#include "const.h"
#include "polevl.h"

namespace xsf {
namespace cephes {

    namespace detail {

        constexpr double ellpk_P[] = {1.37982864606273237150E-4, 2.28025724005875567385E-3, 7.97404013220415179367E-3,
                                      9.85821379021226008714E-3, 6.87489687449949877925E-3, 6.18901033637687613229E-3,
                                      8.79078273952743772254E-3, 1.49380448916805252718E-2, 3.08851465246711995998E-2,
                                      9.65735902811690126535E-2, 1.38629436111989062502E0};

        constexpr double ellpk_Q[] = {2.94078955048598507511E-5, 9.14184723865917226571E-4, 5.94058303753167793257E-3,
                                      1.54850516649762399335E-2, 2.39089602715924892727E-2, 3.01204715227604046988E-2,
                                      3.73774314173823228969E-2, 4.88280347570998239232E-2, 7.03124996963957469739E-2,
                                      1.24999999999870820058E-1, 4.99999999999999999821E-1};

        constexpr double ellpk_C1 = 1.3862943611198906188E0; /* log(4) */

    } // namespace detail

    XSF_HOST_DEVICE inline double ellpk(double x) {

        if (x < 0.0) {
            set_error("ellpk", SF_ERROR_DOMAIN, NULL);
            return (std::numeric_limits<double>::quiet_NaN());
        }

        if (x > 1.0) {
            if (std::isinf(x)) {
                return 0.0;
            }
            return ellpk(1 / x) / std::sqrt(x);
        }

        if (x > detail::MACHEP) {
            return (polevl(x, detail::ellpk_P, 10) - std::log(x) * polevl(x, detail::ellpk_Q, 10));
        } else {
            if (x == 0.0) {
                set_error("ellpk", SF_ERROR_SINGULAR, NULL);
                return (std::numeric_limits<double>::infinity());
            } else {
                return (detail::ellpk_C1 - 0.5 * std::log(x));
            }
        }
    }
} // namespace cephes
} // namespace xsf
