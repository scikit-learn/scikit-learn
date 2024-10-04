/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     k1.c
 *
 *     Modified Bessel function, third kind, order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k1();
 *
 * y = k1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the modified Bessel function of the third kind
 * of order one of the argument.
 *
 * The range is partitioned into the two intervals [0,2] and
 * (2, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       1.2e-15     1.6e-16
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * k1 domain          x <= 0          INFINITY
 *
 */
/*							k1e.c
 *
 *	Modified Bessel function, third kind, order one,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k1e();
 *
 * y = k1e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of the third kind of order one of the argument:
 *
 *      k1e(x) = exp(x) * k1(x).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       7.8e-16     1.2e-16
 * See k1().
 *
 */

/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 2000 by Stephen L. Moshier
 */
#pragma once

#include "../config.h"
#include "../error.h"

#include "chbevl.h"
#include "const.h"

namespace special {
namespace cephes {

    namespace detail {
        /* Chebyshev coefficients for x(K1(x) - log(x/2) I1(x))
         * in the interval [0,2].
         *
         * lim(x->0){ x(K1(x) - log(x/2) I1(x)) } = 1.
         */

        constexpr double k1_A[] = {
            -7.02386347938628759343E-18, -2.42744985051936593393E-15, -6.66690169419932900609E-13,
            -1.41148839263352776110E-10, -2.21338763073472585583E-8,  -2.43340614156596823496E-6,
            -1.73028895751305206302E-4,  -6.97572385963986435018E-3,  -1.22611180822657148235E-1,
            -3.53155960776544875667E-1,  1.52530022733894777053E0};

        /* Chebyshev coefficients for exp(x) sqrt(x) K1(x)
         * in the interval [2,infinity].
         *
         * lim(x->inf){ exp(x) sqrt(x) K1(x) } = sqrt(pi/2).
         */
        constexpr double k1_B[] = {
            -5.75674448366501715755E-18, 1.79405087314755922667E-17,  -5.68946255844285935196E-17,
            1.83809354436663880070E-16,  -6.05704724837331885336E-16, 2.03870316562433424052E-15,
            -7.01983709041831346144E-15, 2.47715442448130437068E-14,  -8.97670518232499435011E-14,
            3.34841966607842919884E-13,  -1.28917396095102890680E-12, 5.13963967348173025100E-12,
            -2.12996783842756842877E-11, 9.21831518760500529508E-11,  -4.19035475934189648750E-10,
            2.01504975519703286596E-9,   -1.03457624656780970260E-8,  5.74108412545004946722E-8,
            -3.50196060308781257119E-7,  2.40648494783721712015E-6,   -1.93619797416608296024E-5,
            1.95215518471351631108E-4,   -2.85781685962277938680E-3,  1.03923736576817238437E-1,
            2.72062619048444266945E0};

    } // namespace detail

    SPECFUN_HOST_DEVICE inline double k1(double x) {
        double y, z;

        if (x == 0.0) {
            set_error("k1", SF_ERROR_SINGULAR, NULL);
            return std::numeric_limits<double>::infinity();
        } else if (x < 0.0) {
            set_error("k1", SF_ERROR_DOMAIN, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }
        z = 0.5 * x;

        if (x <= 2.0) {
            y = x * x - 2.0;
            y = std::log(z) * i1(x) + chbevl(y, detail::k1_A, 11) / x;
            return (y);
        }

        return (std::exp(-x) * chbevl(8.0 / x - 2.0, detail::k1_B, 25) / std::sqrt(x));
    }

    SPECFUN_HOST_DEVICE double k1e(double x) {
        double y;

        if (x == 0.0) {
            set_error("k1e", SF_ERROR_SINGULAR, NULL);
            return std::numeric_limits<double>::infinity();
        } else if (x < 0.0) {
            set_error("k1e", SF_ERROR_DOMAIN, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }

        if (x <= 2.0) {
            y = x * x - 2.0;
            y = std::log(0.5 * x) * i1(x) + chbevl(y, detail::k1_A, 11) / x;
            return (y * exp(x));
        }

        return (chbevl(8.0 / x - 2.0, detail::k1_B, 25) / std::sqrt(x));
    }

} // namespace cephes
} // namespace special
