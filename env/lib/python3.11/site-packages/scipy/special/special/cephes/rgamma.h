/*                                             rgamma.c
 *
 *     Reciprocal Gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, rgamma();
 *
 * y = rgamma( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns one divided by the Gamma function of the argument.
 *
 * The function is approximated by a Chebyshev expansion in
 * the interval [0,1].  Range reduction is by recurrence
 * for arguments between -34.034 and +34.84425627277176174.
 * 0 is returned for positive arguments outside this
 * range.  For arguments less than -34.034 the cosecant
 * reflection formula is applied; lograrithms are employed
 * to avoid unnecessary overflow.
 *
 * The reciprocal Gamma function has no singularities,
 * but overflow and underflow may occur for large arguments.
 * These conditions return either INFINITY or 0 with
 * appropriate sign.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -30,+30      30000       1.1e-15     2.0e-16
 * For arguments less than -34.034 the peak error is on the
 * order of 5e-15 (DEC), excepting overflow or underflow.
 */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1985, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
#pragma once

#include "../config.h"
#include "../error.h"
#include "chbevl.h"
#include "const.h"
#include "gamma.h"
#include "trig.h"

namespace special {
namespace cephes {

    namespace detail {

        /* Chebyshev coefficients for reciprocal Gamma function
         * in interval 0 to 1.  Function is 1/(x Gamma(x)) - 1
         */

        constexpr double rgamma_R[] = {
            3.13173458231230000000E-17, -6.70718606477908000000E-16, 2.20039078172259550000E-15,
            2.47691630348254132600E-13, -6.60074100411295197440E-12, 5.13850186324226978840E-11,
            1.08965386454418662084E-9,  -3.33964630686836942556E-8,  2.68975996440595483619E-7,
            2.96001177518801696639E-6,  -8.04814124978471142852E-5,  4.16609138709688864714E-4,
            5.06579864028608725080E-3,  -6.41925436109158228810E-2,  -4.98558728684003594785E-3,
            1.27546015610523951063E-1};

    } // namespace detail

    SPECFUN_HOST_DEVICE double rgamma(double x) {
        double w, y, z;
        int sign;

        if (x > 34.84425627277176174) {
            return std::exp(-lgam(x));
        }
        if (x < -34.034) {
            w = -x;
            z = sinpi(w);
            if (z == 0.0) {
                return 0.0;
            }
            if (z < 0.0) {
                sign = 1;
                z = -z;
            } else {
                sign = -1;
            }

            y = std::log(w * z) - std::log(M_PI) + lgam(w);
            if (y < -detail::MAXLOG) {
                set_error("rgamma", SF_ERROR_UNDERFLOW, NULL);
                return (sign * 0.0);
            }
            if (y > detail::MAXLOG) {
                set_error("rgamma", SF_ERROR_OVERFLOW, NULL);
                return (sign * std::numeric_limits<double>::infinity());
            }
            return (sign * std::exp(y));
        }
        z = 1.0;
        w = x;

        while (w > 1.0) { /* Downward recurrence */
            w -= 1.0;
            z *= w;
        }
        while (w < 0.0) { /* Upward recurrence */
            z /= w;
            w += 1.0;
        }
        if (w == 0.0) /* Nonpositive integer */
            return (0.0);
        if (w == 1.0) /* Other integer */
            return (1.0 / z);

        y = w * (1.0 + chbevl(4.0 * w - 2.0, detail::rgamma_R, 16)) / z;
        return (y);
    }
} // namespace cephes
} // namespace special
