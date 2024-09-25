/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     i1.c
 *
 *     Modified Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i1();
 *
 * y = i1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order one of the
 * argument.
 *
 * The function is defined as i1(x) = -i j1( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       1.9e-15     2.1e-16
 *
 *
 */
/*							i1e.c
 *
 *	Modified Bessel function of order one,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i1e();
 *
 * y = i1e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of order one of the argument.
 *
 * The function is defined as i1(x) = -i exp(-|x|) j1( ix ).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       2.0e-15     2.0e-16
 * See i1().
 *
 */

/*                                                     i1.c 2          */

/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1985, 1987, 2000 by Stephen L. Moshier
 */
#pragma once

#include "../config.h"
#include "chbevl.h"

namespace special {
namespace cephes {

    namespace detail {

        /* Chebyshev coefficients for exp(-x) I1(x) / x
         * in the interval [0,8].
         *
         * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
         */

        constexpr double i1_A[] = {
            2.77791411276104639959E-18,  -2.11142121435816608115E-17, 1.55363195773620046921E-16,
            -1.10559694773538630805E-15, 7.60068429473540693410E-15,  -5.04218550472791168711E-14,
            3.22379336594557470981E-13,  -1.98397439776494371520E-12, 1.17361862988909016308E-11,
            -6.66348972350202774223E-11, 3.62559028155211703701E-10,  -1.88724975172282928790E-9,
            9.38153738649577178388E-9,   -4.44505912879632808065E-8,  2.00329475355213526229E-7,
            -8.56872026469545474066E-7,  3.47025130813767847674E-6,   -1.32731636560394358279E-5,
            4.78156510755005422638E-5,   -1.61760815825896745588E-4,  5.12285956168575772895E-4,
            -1.51357245063125314899E-3,  4.15642294431288815669E-3,   -1.05640848946261981558E-2,
            2.47264490306265168283E-2,   -5.29459812080949914269E-2,  1.02643658689847095384E-1,
            -1.76416518357834055153E-1,  2.52587186443633654823E-1};

        /* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
         * in the inverted interval [8,infinity].
         *
         * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
         */
        constexpr double i1_B[] = {
            7.51729631084210481353E-18,  4.41434832307170791151E-18,  -4.65030536848935832153E-17,
            -3.20952592199342395980E-17, 2.96262899764595013876E-16,  3.30820231092092828324E-16,
            -1.88035477551078244854E-15, -3.81440307243700780478E-15, 1.04202769841288027642E-14,
            4.27244001671195135429E-14,  -2.10154184277266431302E-14, -4.08355111109219731823E-13,
            -7.19855177624590851209E-13, 2.03562854414708950722E-12,  1.41258074366137813316E-11,
            3.25260358301548823856E-11,  -1.89749581235054123450E-11, -5.58974346219658380687E-10,
            -3.83538038596423702205E-9,  -2.63146884688951950684E-8,  -2.51223623787020892529E-7,
            -3.88256480887769039346E-6,  -1.10588938762623716291E-4,  -9.76109749136146840777E-3,
            7.78576235018280120474E-1};

    } // namespace detail

    SPECFUN_HOST_DEVICE inline double i1(double x) {
        double y, z;

        z = std::abs(x);
        if (z <= 8.0) {
            y = (z / 2.0) - 2.0;
            z = chbevl(y, detail::i1_A, 29) * z * std::exp(z);
        } else {
            z = std::exp(z) * chbevl(32.0 / z - 2.0, detail::i1_B, 25) / std::sqrt(z);
        }
        if (x < 0.0)
            z = -z;
        return (z);
    }

    /*                                                     i1e()   */

    SPECFUN_HOST_DEVICE inline double i1e(double x) {
        double y, z;

        z = std::abs(x);
        if (z <= 8.0) {
            y = (z / 2.0) - 2.0;
            z = chbevl(y, detail::i1_A, 29) * z;
        } else {
            z = chbevl(32.0 / z - 2.0, detail::i1_B, 25) / std::sqrt(z);
        }
        if (x < 0.0)
            z = -z;
        return (z);
    }

} // namespace cephes
} // namespace special
