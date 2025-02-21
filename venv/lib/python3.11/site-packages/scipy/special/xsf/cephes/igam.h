/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     igam.c
 *
 *     Incomplete Gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, y, igam();
 *
 * y = igam( a, x );
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *                           x
 *                            -
 *                   1       | |  -t  a-1
 *  igam(a,x)  =   -----     |   e   t   dt.
 *                  -      | |
 *                 | (a)    -
 *                           0
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30       200000       3.6e-14     2.9e-15
 *    IEEE      0,100      300000       9.9e-14     1.5e-14
 */
/*							igamc()
 *
 *	Complemented incomplete Gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, y, igamc();
 *
 * y = igamc( a, x );
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *
 *  igamc(a,x)   =   1 - igam(a,x)
 *
 *                            inf.
 *                              -
 *                     1       | |  -t  a-1
 *               =   -----     |   e   t   dt.
 *                    -      | |
 *                   | (a)    -
 *                             x
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 * ACCURACY:
 *
 * Tested at random a, x.
 *                a         x                      Relative error:
 * arithmetic   domain   domain     # trials      peak         rms
 *    IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
 *    IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15
 */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1985, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

/* Sources
 * [1] "The Digital Library of Mathematical Functions", dlmf.nist.gov
 * [2] Maddock et. al., "Incomplete Gamma Functions",
 *     https://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/sf_gamma/igamma.html
 */

/* Scipy changes:
 * - 05-01-2016: added asymptotic expansion for igam to improve the
 *   a ~ x regime.
 * - 06-19-2016: additional series expansion added for igamc to
 *   improve accuracy at small arguments.
 * - 06-24-2016: better choice of domain for the asymptotic series;
 *   improvements in accuracy for the asymptotic series when a and x
 *   are very close.
 */
#pragma once

#include "../config.h"
#include "../error.h"

#include "const.h"
#include "gamma.h"
#include "igam_asymp_coeff.h"
#include "lanczos.h"
#include "ndtr.h"
#include "unity.h"

namespace xsf {
namespace cephes {

    namespace detail {

        constexpr int igam_MAXITER = 2000;
        constexpr int IGAM = 1;
        constexpr int IGAMC = 0;
        constexpr double igam_SMALL = 20;
        constexpr double igam_LARGE = 200;
        constexpr double igam_SMALLRATIO = 0.3;
        constexpr double igam_LARGERATIO = 4.5;

        constexpr double igam_big = 4.503599627370496e15;
        constexpr double igam_biginv = 2.22044604925031308085e-16;

        /* Compute
         *
         * x^a * exp(-x) / gamma(a)
         *
         * corrected from (15) and (16) in [2] by replacing exp(x - a) with
         * exp(a - x).
         */
        XSF_HOST_DEVICE inline double igam_fac(double a, double x) {
            double ax, fac, res, num;

            if (std::abs(a - x) > 0.4 * std::abs(a)) {
                ax = a * std::log(x) - x - xsf::cephes::lgam(a);
                if (ax < -MAXLOG) {
                    set_error("igam", SF_ERROR_UNDERFLOW, NULL);
                    return 0.0;
                }
                return std::exp(ax);
            }

            fac = a + xsf::cephes::lanczos_g - 0.5;
            res = std::sqrt(fac / std::exp(1)) / xsf::cephes::lanczos_sum_expg_scaled(a);

            if ((a < 200) && (x < 200)) {
                res *= std::exp(a - x) * std::pow(x / fac, a);
            } else {
                num = x - a - xsf::cephes::lanczos_g + 0.5;
                res *= std::exp(a * xsf::cephes::log1pmx(num / fac) + x * (0.5 - xsf::cephes::lanczos_g) / fac);
            }

            return res;
        }

        /* Compute igamc using DLMF 8.9.2. */
        XSF_HOST_DEVICE inline double igamc_continued_fraction(double a, double x) {
            int i;
            double ans, ax, c, yc, r, t, y, z;
            double pk, pkm1, pkm2, qk, qkm1, qkm2;

            ax = igam_fac(a, x);
            if (ax == 0.0) {
                return 0.0;
            }

            /* continued fraction */
            y = 1.0 - a;
            z = x + y + 1.0;
            c = 0.0;
            pkm2 = 1.0;
            qkm2 = x;
            pkm1 = x + 1.0;
            qkm1 = z * x;
            ans = pkm1 / qkm1;

            for (i = 0; i < igam_MAXITER; i++) {
                c += 1.0;
                y += 1.0;
                z += 2.0;
                yc = y * c;
                pk = pkm1 * z - pkm2 * yc;
                qk = qkm1 * z - qkm2 * yc;
                if (qk != 0) {
                    r = pk / qk;
                    t = std::abs((ans - r) / r);
                    ans = r;
                } else
                    t = 1.0;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                if (std::abs(pk) > igam_big) {
                    pkm2 *= igam_biginv;
                    pkm1 *= igam_biginv;
                    qkm2 *= igam_biginv;
                    qkm1 *= igam_biginv;
                }
                if (t <= MACHEP) {
                    break;
                }
            }

            return (ans * ax);
        }

        /* Compute igam using DLMF 8.11.4. */
        XSF_HOST_DEVICE inline double igam_series(double a, double x) {
            int i;
            double ans, ax, c, r;

            ax = igam_fac(a, x);
            if (ax == 0.0) {
                return 0.0;
            }

            /* power series */
            r = a;
            c = 1.0;
            ans = 1.0;

            for (i = 0; i < igam_MAXITER; i++) {
                r += 1.0;
                c *= x / r;
                ans += c;
                if (c <= MACHEP * ans) {
                    break;
                }
            }

            return (ans * ax / a);
        }

        /* Compute igamc using DLMF 8.7.3. This is related to the series in
         * igam_series but extra care is taken to avoid cancellation.
         */
        XSF_HOST_DEVICE inline double igamc_series(double a, double x) {
            int n;
            double fac = 1;
            double sum = 0;
            double term, logx;

            for (n = 1; n < igam_MAXITER; n++) {
                fac *= -x / n;
                term = fac / (a + n);
                sum += term;
                if (std::abs(term) <= MACHEP * std::abs(sum)) {
                    break;
                }
            }

            logx = std::log(x);
            term = -xsf::cephes::expm1(a * logx - xsf::cephes::lgam1p(a));
            return term - std::exp(a * logx - xsf::cephes::lgam(a)) * sum;
        }

        /* Compute igam/igamc using DLMF 8.12.3/8.12.4. */
        XSF_HOST_DEVICE inline double asymptotic_series(double a, double x, int func) {
            int k, n, sgn;
            int maxpow = 0;
            double lambda = x / a;
            double sigma = (x - a) / a;
            double eta, res, ck, ckterm, term, absterm;
            double absoldterm = std::numeric_limits<double>::infinity();
            double etapow[detail::igam_asymp_coeff_N] = {1};
            double sum = 0;
            double afac = 1;

            if (func == detail::IGAM) {
                sgn = -1;
            } else {
                sgn = 1;
            }

            if (lambda > 1) {
                eta = std::sqrt(-2 * xsf::cephes::log1pmx(sigma));
            } else if (lambda < 1) {
                eta = -std::sqrt(-2 * xsf::cephes::log1pmx(sigma));
            } else {
                eta = 0;
            }
            res = 0.5 * xsf::cephes::erfc(sgn * eta * std::sqrt(a / 2));

            for (k = 0; k < igam_asymp_coeff_K; k++) {
                ck = igam_asymp_coeff_d[k][0];
                for (n = 1; n < igam_asymp_coeff_N; n++) {
                    if (n > maxpow) {
                        etapow[n] = eta * etapow[n - 1];
                        maxpow += 1;
                    }
                    ckterm = igam_asymp_coeff_d[k][n] * etapow[n];
                    ck += ckterm;
                    if (std::abs(ckterm) < MACHEP * std::abs(ck)) {
                        break;
                    }
                }
                term = ck * afac;
                absterm = std::abs(term);
                if (absterm > absoldterm) {
                    break;
                }
                sum += term;
                if (absterm < MACHEP * std::abs(sum)) {
                    break;
                }
                absoldterm = absterm;
                afac /= a;
            }
            res += sgn * std::exp(-0.5 * a * eta * eta) * sum / std::sqrt(2 * M_PI * a);

            return res;
        }

    } // namespace detail

    XSF_HOST_DEVICE inline double igamc(double a, double x);

    XSF_HOST_DEVICE inline double igam(double a, double x) {
        double absxma_a;

        if (x < 0 || a < 0) {
            set_error("gammainc", SF_ERROR_DOMAIN, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        } else if (a == 0) {
            if (x > 0) {
                return 1;
            } else {
                return std::numeric_limits<double>::quiet_NaN();
            }
        } else if (x == 0) {
            /* Zero integration limit */
            return 0;
        } else if (std::isinf(a)) {
            if (std::isinf(x)) {
                return std::numeric_limits<double>::quiet_NaN();
            }
            return 0;
        } else if (std::isinf(x)) {
            return 1;
        }

        /* Asymptotic regime where a ~ x; see [2]. */
        absxma_a = std::abs(x - a) / a;
        if ((a > detail::igam_SMALL) && (a < detail::igam_LARGE) && (absxma_a < detail::igam_SMALLRATIO)) {
            return detail::asymptotic_series(a, x, detail::IGAM);
        } else if ((a > detail::igam_LARGE) && (absxma_a < detail::igam_LARGERATIO / std::sqrt(a))) {
            return detail::asymptotic_series(a, x, detail::IGAM);
        }

        if ((x > 1.0) && (x > a)) {
            return (1.0 - igamc(a, x));
        }

        return detail::igam_series(a, x);
    }

    XSF_HOST_DEVICE double igamc(double a, double x) {
        double absxma_a;

        if (x < 0 || a < 0) {
            set_error("gammaincc", SF_ERROR_DOMAIN, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        } else if (a == 0) {
            if (x > 0) {
                return 0;
            } else {
                return std::numeric_limits<double>::quiet_NaN();
            }
        } else if (x == 0) {
            return 1;
        } else if (std::isinf(a)) {
            if (std::isinf(x)) {
                return std::numeric_limits<double>::quiet_NaN();
            }
            return 1;
        } else if (std::isinf(x)) {
            return 0;
        }

        /* Asymptotic regime where a ~ x; see [2]. */
        absxma_a = std::abs(x - a) / a;
        if ((a > detail::igam_SMALL) && (a < detail::igam_LARGE) && (absxma_a < detail::igam_SMALLRATIO)) {
            return detail::asymptotic_series(a, x, detail::IGAMC);
        } else if ((a > detail::igam_LARGE) && (absxma_a < detail::igam_LARGERATIO / std::sqrt(a))) {
            return detail::asymptotic_series(a, x, detail::IGAMC);
        }

        /* Everywhere else; see [2]. */
        if (x > 1.1) {
            if (x < a) {
                return 1.0 - detail::igam_series(a, x);
            } else {
                return detail::igamc_continued_fraction(a, x);
            }
        } else if (x <= 0.5) {
            if (-0.4 / std::log(x) < a) {
                return 1.0 - detail::igam_series(a, x);
            } else {
                return detail::igamc_series(a, x);
            }
        } else {
            if (x * 1.1 < a) {
                return 1.0 - detail::igam_series(a, x);
            } else {
                return detail::igamc_series(a, x);
            }
        }
    }

} // namespace cephes
} // namespace xsf
