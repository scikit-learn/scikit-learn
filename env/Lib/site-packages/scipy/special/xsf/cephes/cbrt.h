/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     cbrt.c
 *
 *     Cube root
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, cbrt();
 *
 * y = cbrt( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the cube root of the argument, which may be negative.
 *
 * Range reduction involves determining the power of 2 of
 * the argument.  A polynomial of degree 2 applied to the
 * mantissa, and multiplication by the cube root of 1, 2, or 4
 * approximates the root to within about 0.1%.  Then Newton's
 * iteration is used three times to converge to an accurate
 * result.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,1e308     30000      1.5e-16     5.0e-17
 *
 */
/*							cbrt.c  */

/*
 * Cephes Math Library Release 2.2:  January, 1991
 * Copyright 1984, 1991 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
#pragma once

#include "../config.h"

namespace xsf {
namespace cephes {

    namespace detail {

        constexpr double CBRT2 = 1.2599210498948731647672;
        constexpr double CBRT4 = 1.5874010519681994747517;
        constexpr double CBRT2I = 0.79370052598409973737585;
        constexpr double CBRT4I = 0.62996052494743658238361;

        XSF_HOST_DEVICE inline double cbrt(double x) {
            int e, rem, sign;
            double z;

            if (!std::isfinite(x)) {
                return x;
            }
            if (x == 0) {
                return (x);
            }
            if (x > 0) {
                sign = 1;
            } else {
                sign = -1;
                x = -x;
            }

            z = x;
            /* extract power of 2, leaving
             * mantissa between 0.5 and 1
             */
            x = std::frexp(x, &e);

            /* Approximate cube root of number between .5 and 1,
             * peak relative error = 9.2e-6
             */
            x = (((-1.3466110473359520655053e-1 * x + 5.4664601366395524503440e-1) * x - 9.5438224771509446525043e-1) *
                     x +
                 1.1399983354717293273738e0) *
                    x +
                4.0238979564544752126924e-1;

            /* exponent divided by 3 */
            if (e >= 0) {
                rem = e;
                e /= 3;
                rem -= 3 * e;
                if (rem == 1) {
                    x *= CBRT2;
                } else if (rem == 2) {
                    x *= CBRT4;
                }
            }
            /* argument less than 1 */
            else {
                e = -e;
                rem = e;
                e /= 3;
                rem -= 3 * e;
                if (rem == 1) {
                    x *= CBRT2I;
                } else if (rem == 2) {
                    x *= CBRT4I;
                }
                e = -e;
            }

            /* multiply by power of 2 */
            x = std::ldexp(x, e);

            /* Newton iteration */
            x -= (x - (z / (x * x))) * 0.33333333333333333333;
            x -= (x - (z / (x * x))) * 0.33333333333333333333;

            if (sign < 0)
                x = -x;
            return (x);
        }
    } // namespace detail

} // namespace cephes
} // namespace xsf
