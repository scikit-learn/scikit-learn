/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     expn.c
 *
 *             Exponential integral En
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * double x, y, expn();
 *
 * y = expn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the exponential integral
 *
 *                 inf.
 *                   -
 *                  | |   -xt
 *                  |    e
 *      E (x)  =    |    ----  dt.
 *       n          |      n
 *                | |     t
 *                 -
 *                  1
 *
 *
 * Both n and x must be nonnegative.
 *
 * The routine employs either a power series, a continued
 * fraction, or an asymptotic formula depending on the
 * relative values of n and x.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       10000       1.7e-15     3.6e-16
 *
 */

/*                                                     expn.c  */

/* Cephes Math Library Release 1.1:  March, 1985
 * Copyright 1985 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140 */

/* Sources
 * [1] NIST, "The Digital Library of Mathematical Functions", dlmf.nist.gov
 */

/* Scipy changes:
 * - 09-10-2016: improved asymptotic expansion for large n
 */

#pragma once

#include "../config.h"
#include "../error.h"
#include "const.h"
#include "rgamma.h"
#include "polevl.h"

namespace xsf {
namespace cephes {

    namespace detail {

        constexpr int expn_nA = 13;
        constexpr double expn_A0[] = {1.00000000000000000};
        constexpr double expn_A1[] = {1.00000000000000000};
        constexpr double expn_A2[] = {-2.00000000000000000, 1.00000000000000000};
        constexpr double expn_A3[] = {6.00000000000000000, -8.00000000000000000, 1.00000000000000000};
        constexpr double expn_A4[] = {-24.0000000000000000, 58.0000000000000000, -22.0000000000000000,
                                      1.00000000000000000};
        constexpr double expn_A5[] = {120.000000000000000, -444.000000000000000, 328.000000000000000,
                                      -52.0000000000000000, 1.00000000000000000};
        constexpr double expn_A6[] = {-720.000000000000000, 3708.00000000000000,  -4400.00000000000000,
                                      1452.00000000000000,  -114.000000000000000, 1.00000000000000000};
        constexpr double expn_A7[] = {5040.00000000000000,  -33984.0000000000000, 58140.0000000000000,
                                      -32120.0000000000000, 5610.00000000000000,  -240.000000000000000,
                                      1.00000000000000000};
        constexpr double expn_A8[] = {-40320.0000000000000, 341136.000000000000,  -785304.000000000000,
                                      644020.000000000000,  -195800.000000000000, 19950.0000000000000,
                                      -494.000000000000000, 1.00000000000000000};
        constexpr double expn_A9[] = {362880.000000000000,  -3733920.00000000000, 11026296.0000000000,
                                      -12440064.0000000000, 5765500.00000000000,  -1062500.00000000000,
                                      67260.0000000000000,  -1004.00000000000000, 1.00000000000000000};
        constexpr double expn_A10[] = {-3628800.00000000000, 44339040.0000000000,  -162186912.000000000,
                                       238904904.000000000,  -155357384.000000000, 44765000.0000000000,
                                       -5326160.00000000000, 218848.000000000000,  -2026.00000000000000,
                                       1.00000000000000000};
        constexpr double expn_A11[] = {39916800.0000000000,  -568356480.000000000, 2507481216.00000000,
                                       -4642163952.00000000, 4002695088.00000000,  -1648384304.00000000,
                                       314369720.000000000,  -25243904.0000000000, 695038.000000000000,
                                       -4072.00000000000000, 1.00000000000000000};
        constexpr double expn_A12[] = {-479001600.000000000, 7827719040.00000000,  -40788301824.0000000,
                                       92199790224.0000000,  -101180433024.000000, 56041398784.0000000,
                                       -15548960784.0000000, 2051482776.00000000,  -114876376.000000000,
                                       2170626.00000000000,  -8166.00000000000000, 1.00000000000000000};
        constexpr const double *expn_A[] = {expn_A0, expn_A1, expn_A2, expn_A3,  expn_A4,  expn_A5, expn_A6,
                                            expn_A7, expn_A8, expn_A9, expn_A10, expn_A11, expn_A12};
        constexpr int expn_Adegs[] = {0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

        /* Asymptotic expansion for large n, DLMF 8.20(ii) */
        XSF_HOST_DEVICE double expn_large_n(int n, double x) {
            int k;
            double p = n;
            double lambda = x / p;
            double multiplier = 1 / p / (lambda + 1) / (lambda + 1);
            double fac = 1;
            double res = 1; /* A[0] = 1 */
            double expfac, term;

            expfac = std::exp(-lambda * p) / (lambda + 1) / p;
            if (expfac == 0) {
                set_error("expn", SF_ERROR_UNDERFLOW, NULL);
                return 0;
            }

            /* Do the k = 1 term outside the loop since A[1] = 1 */
            fac *= multiplier;
            res += fac;

            for (k = 2; k < expn_nA; k++) {
                fac *= multiplier;
                term = fac * polevl(lambda, expn_A[k], expn_Adegs[k]);
                res += term;
                if (std::abs(term) < MACHEP * std::abs(res)) {
                    break;
                }
            }

            return expfac * res;
        }
    } // namespace detail

    XSF_HOST_DEVICE double expn(int n, double x) {
        double ans, r, t, yk, xk;
        double pk, pkm1, pkm2, qk, qkm1, qkm2;
        double psi, z;
        int i, k;
        constexpr double big = 1.44115188075855872E+17;

        if (std::isnan(x)) {
            return std::numeric_limits<double>::quiet_NaN();
        } else if (n < 0 || x < 0) {
            set_error("expn", SF_ERROR_DOMAIN, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }

        if (x > detail::MAXLOG) {
            return (0.0);
        }

        if (x == 0.0) {
            if (n < 2) {
                set_error("expn", SF_ERROR_SINGULAR, NULL);
                return std::numeric_limits<double>::infinity();
            } else {
                return (1.0 / (n - 1.0));
            }
        }

        if (n == 0) {
            return (std::exp(-x) / x);
        }

        /* Asymptotic expansion for large n, DLMF 8.20(ii) */
        if (n > 50) {
            ans = detail::expn_large_n(n, x);
            return ans;
        }

        /* Continued fraction, DLMF 8.19.17 */
        if (x > 1.0) {
            k = 1;
            pkm2 = 1.0;
            qkm2 = x;
            pkm1 = 1.0;
            qkm1 = x + n;
            ans = pkm1 / qkm1;

            do {
                k += 1;
                if (k & 1) {
                    yk = 1.0;
                    xk = n + (k - 1) / 2;
                } else {
                    yk = x;
                    xk = k / 2;
                }
                pk = pkm1 * yk + pkm2 * xk;
                qk = qkm1 * yk + qkm2 * xk;
                if (qk != 0) {
                    r = pk / qk;
                    t = std::abs((ans - r) / r);
                    ans = r;
                } else {
                    t = 1.0;
                }
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                if (std::abs(pk) > big) {
                    pkm2 /= big;
                    pkm1 /= big;
                    qkm2 /= big;
                    qkm1 /= big;
                }
            } while (t > detail::MACHEP);

            ans *= std::exp(-x);
            return ans;
        }

        /* Power series expansion, DLMF 8.19.8 */
        psi = -detail::SCIPY_EULER - std::log(x);
        for (i = 1; i < n; i++) {
            psi = psi + 1.0 / i;
        }

        z = -x;
        xk = 0.0;
        yk = 1.0;
        pk = 1.0 - n;
        if (n == 1) {
            ans = 0.0;
        } else {
            ans = 1.0 / pk;
        }
        do {
            xk += 1.0;
            yk *= z / xk;
            pk += 1.0;
            if (pk != 0.0) {
                ans += yk / pk;
            }
            if (ans != 0.0)
                t = std::abs(yk / ans);
            else
                t = 1.0;
        } while (t > detail::MACHEP);
        k = xk;
        t = n;
        r = n - 1;
        ans = (std::pow(z, r) * psi * rgamma(t)) - ans;
        return ans;
    }

} // namespace cephes
} // namespace xsf
