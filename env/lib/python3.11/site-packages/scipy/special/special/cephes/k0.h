/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     k0.c
 *
 *     Modified Bessel function, third kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k0();
 *
 * y = k0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of the third kind
 * of order zero of the argument.
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 * Tested at 2000 random points between 0 and 8.  Peak absolute
 * error (relative when K0 > 1) was 1.46e-14; rms, 4.26e-15.
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       1.2e-15     1.6e-16
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 *  K0 domain          x <= 0          INFINITY
 *
 */
/*							k0e()
 *
 *	Modified Bessel function, third kind, order zero,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k0e();
 *
 * y = k0e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of the third kind of order zero of the argument.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       1.4e-15     1.4e-16
 * See k0().
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
#include "i0.h"

namespace special {
namespace cephes {

    namespace detail {
        /* Chebyshev coefficients for K0(x) + log(x/2) I0(x)
         * in the interval [0,2].  The odd order coefficients are all
         * zero; only the even order coefficients are listed.
         *
         * lim(x->0){ K0(x) + log(x/2) I0(x) } = -EUL.
         */

        constexpr double k0_A[] = {1.37446543561352307156E-16, 4.25981614279661018399E-14, 1.03496952576338420167E-11,
                                   1.90451637722020886025E-9,  2.53479107902614945675E-7,  2.28621210311945178607E-5,
                                   1.26461541144692592338E-3,  3.59799365153615016266E-2,  3.44289899924628486886E-1,
                                   -5.35327393233902768720E-1};

        /* Chebyshev coefficients for exp(x) sqrt(x) K0(x)
         * in the inverted interval [2,infinity].
         *
         * lim(x->inf){ exp(x) sqrt(x) K0(x) } = sqrt(pi/2).
         */
        constexpr double k0_B[] = {
            5.30043377268626276149E-18,  -1.64758043015242134646E-17, 5.21039150503902756861E-17,
            -1.67823109680541210385E-16, 5.51205597852431940784E-16,  -1.84859337734377901440E-15,
            6.34007647740507060557E-15,  -2.22751332699166985548E-14, 8.03289077536357521100E-14,
            -2.98009692317273043925E-13, 1.14034058820847496303E-12,  -4.51459788337394416547E-12,
            1.85594911495471785253E-11,  -7.95748924447710747776E-11, 3.57739728140030116597E-10,
            -1.69753450938905987466E-9,  8.57403401741422608519E-9,   -4.66048989768794782956E-8,
            2.76681363944501510342E-7,   -1.83175552271911948767E-6,  1.39498137188764993662E-5,
            -1.28495495816278026384E-4,  1.56988388573005337491E-3,   -3.14481013119645005427E-2,
            2.44030308206595545468E0};

    } // namespace detail

    SPECFUN_HOST_DEVICE inline double k0(double x) {
        double y, z;

        if (x == 0.0) {
            set_error("k0", SF_ERROR_SINGULAR, NULL);
            return std::numeric_limits<double>::infinity();
        } else if (x < 0.0) {
            set_error("k0", SF_ERROR_DOMAIN, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }

        if (x <= 2.0) {
            y = x * x - 2.0;
            y = chbevl(y, detail::k0_A, 10) - std::log(0.5 * x) * i0(x);
            return (y);
        }
        z = 8.0 / x - 2.0;
        y = std::exp(-x) * chbevl(z, detail::k0_B, 25) / std::sqrt(x);
        return (y);
    }

    SPECFUN_HOST_DEVICE double inline k0e(double x) {
        double y;

        if (x == 0.0) {
            set_error("k0e", SF_ERROR_SINGULAR, NULL);
            return std::numeric_limits<double>::infinity();
        } else if (x < 0.0) {
            set_error("k0e", SF_ERROR_DOMAIN, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }

        if (x <= 2.0) {
            y = x * x - 2.0;
            y = chbevl(y, detail::k0_A, 10) - std::log(0.5 * x) * i0(x);
            return (y * exp(x));
        }

        y = chbevl(8.0 / x - 2.0, detail::k0_B, 25) / std::sqrt(x);
        return (y);
    }

} // namespace cephes
} // namespace special
