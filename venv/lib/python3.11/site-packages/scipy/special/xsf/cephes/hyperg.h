/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     hyperg.c
 *
 *     Confluent hypergeometric function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, hyperg();
 *
 * y = hyperg( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the confluent hypergeometric function
 *
 *                          1           2
 *                       a x    a(a+1) x
 *   F ( a,b;x )  =  1 + ---- + --------- + ...
 *  1 1                  b 1!   b(b+1) 2!
 *
 * Many higher transcendental functions are special cases of
 * this power series.
 *
 * As is evident from the formula, b must not be a negative
 * integer or zero unless a is an integer with 0 >= a > b.
 *
 * The routine attempts both a direct summation of the series
 * and an asymptotic expansion.  In each case error due to
 * roundoff, cancellation, and nonconvergence is estimated.
 * The result with smaller estimated error is returned.
 *
 *
 *
 * ACCURACY:
 *
 * Tested at random points (a, b, x), all three variables
 * ranging from 0 to 30.
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30        30000       1.8e-14     1.1e-15
 *
 * Larger errors can be observed when b is near a negative
 * integer or zero.  Certain combinations of arguments yield
 * serious cancellation error in the power series summation
 * and also are not in the region of near convergence of the
 * asymptotic series.  An error message is printed if the
 * self-estimated relative error is greater than 1.0e-12.
 *
 */

/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 1988, 2000 by Stephen L. Moshier
 */

#pragma once

#include "../config.h"
#include "../error.h"

#include "const.h"
#include "gamma.h"
#include "rgamma.h"

namespace xsf {
namespace cephes {

    namespace detail {

        /* the `type` parameter determines what converging factor to use */
        XSF_HOST_DEVICE inline double hyp2f0(double a, double b, double x, int type, double *err) {
            double a0, alast, t, tlast, maxt;
            double n, an, bn, u, sum, temp;

            an = a;
            bn = b;
            a0 = 1.0e0;
            alast = 1.0e0;
            sum = 0.0;
            n = 1.0e0;
            t = 1.0e0;
            tlast = 1.0e9;
            maxt = 0.0;

            do {
                if (an == 0)
                    goto pdone;
                if (bn == 0)
                    goto pdone;

                u = an * (bn * x / n);

                /* check for blowup */
                temp = std::abs(u);
                if ((temp > 1.0) && (maxt > (std::numeric_limits<double>::max() / temp)))
                    goto error;

                a0 *= u;
                t = std::abs(a0);

                /* terminating condition for asymptotic series:
                 * the series is divergent (if a or b is not a negative integer),
                 * but its leading part can be used as an asymptotic expansion
                 */
                if (t > tlast)
                    goto ndone;

                tlast = t;
                sum += alast; /* the sum is one term behind */
                alast = a0;

                if (n > 200)
                    goto ndone;

                an += 1.0e0;
                bn += 1.0e0;
                n += 1.0e0;
                if (t > maxt)
                    maxt = t;
            } while (t > MACHEP);

        pdone: /* series converged! */

            /* estimate error due to roundoff and cancellation */
            *err = std::abs(MACHEP * (n + maxt));

            alast = a0;
            goto done;

        ndone: /* series did not converge */

            /* The following "Converging factors" are supposed to improve accuracy,
             * but do not actually seem to accomplish very much. */

            n -= 1.0;
            x = 1.0 / x;

            switch (type) { /* "type" given as subroutine argument */
            case 1:
                alast *= (0.5 + (0.125 + 0.25 * b - 0.5 * a + 0.25 * x - 0.25 * n) / x);
                break;

            case 2:
                alast *= 2.0 / 3.0 - b + 2.0 * a + x - n;
                break;

            default:;
            }

            /* estimate error due to roundoff, cancellation, and nonconvergence */
            *err = MACHEP * (n + maxt) + std::abs(a0);

        done:
            sum += alast;
            return (sum);

            /* series blew up: */
        error:
            *err = std::numeric_limits<double>::infinity();
            set_error("hyperg", SF_ERROR_NO_RESULT, NULL);
            return (sum);
        }

        /* asymptotic formula for hypergeometric function:
         *
         *        (    -a
         *  --    ( |z|
         * |  (b) ( -------- 2f0( a, 1+a-b, -1/x )
         *        (  --
         *        ( |  (b-a)
         *
         *
         *                                x    a-b                     )
         *                               e  |x|                        )
         *                             + -------- 2f0( b-a, 1-a, 1/x ) )
         *                                --                           )
         *                               |  (a)                        )
         */

        XSF_HOST_DEVICE inline double hy1f1a(double a, double b, double x, double *err) {
            double h1, h2, t, u, temp, acanc, asum, err1, err2;

            if (x == 0) {
                acanc = 1.0;
                asum = std::numeric_limits<double>::infinity();
                goto adone;
            }
            temp = std::log(std::abs(x));
            t = x + temp * (a - b);
            u = -temp * a;

            if (b > 0) {
                temp = xsf::cephes::lgam(b);
                t += temp;
                u += temp;
            }

            h1 = hyp2f0(a, a - b + 1, -1.0 / x, 1, &err1);

            temp = std::exp(u) * xsf::cephes::rgamma(b - a);
            h1 *= temp;
            err1 *= temp;

            h2 = hyp2f0(b - a, 1.0 - a, 1.0 / x, 2, &err2);

            if (a < 0)
                temp = std::exp(t) * xsf::cephes::rgamma(a);
            else
                temp = std::exp(t - xsf::cephes::lgam(a));

            h2 *= temp;
            err2 *= temp;

            if (x < 0.0)
                asum = h1;
            else
                asum = h2;

            acanc = std::abs(err1) + std::abs(err2);

            if (b < 0) {
                temp = xsf::cephes::Gamma(b);
                asum *= temp;
                acanc *= std::abs(temp);
            }

            if (asum != 0.0)
                acanc /= std::abs(asum);

            if (acanc != acanc)
                /* nan */
                acanc = 1.0;

            if (std::isinf(asum))
                /* infinity */
                acanc = 0;

            acanc *= 30.0; /* fudge factor, since error of asymptotic formula
                            * often seems this much larger than advertised */
        adone:
            *err = acanc;
            return (asum);
        }

        /* Power series summation for confluent hypergeometric function */
        XSF_HOST_DEVICE inline double hy1f1p(double a, double b, double x, double *err) {
            double n, a0, sum, t, u, temp, maxn;
            double an, bn, maxt;
            double y, c, sumc;

            /* set up for power series summation */
            an = a;
            bn = b;
            a0 = 1.0;
            sum = 1.0;
            c = 0.0;
            n = 1.0;
            t = 1.0;
            maxt = 0.0;
            *err = 1.0;

            maxn = 200.0 + 2 * fabs(a) + 2 * fabs(b);

            while (t > MACHEP) {
                if (bn == 0) { /* check bn first since if both   */
                    sf_error("hyperg", SF_ERROR_SINGULAR, NULL);
                    return (std::numeric_limits<double>::infinity()); /* an and bn are zero it is     */
                }
                if (an == 0) /* a singularity            */
                    return (sum);
                if (n > maxn) {
                    /* too many terms; take the last one as error estimate */
                    c = std::abs(c) + std::abs(t) * 50.0;
                    goto pdone;
                }
                u = x * (an / (bn * n));

                /* check for blowup */
                temp = std::abs(u);
                if ((temp > 1.0) && (maxt > (std::numeric_limits<double>::max() / temp))) {
                    *err = 1.0; /* blowup: estimate 100% error */
                    return sum;
                }

                a0 *= u;

                y = a0 - c;
                sumc = sum + y;
                c = (sumc - sum) - y;
                sum = sumc;

                t = std::abs(a0);

                an += 1.0;
                bn += 1.0;
                n += 1.0;
            }

        pdone:

            /* estimate error due to roundoff and cancellation */
            if (sum != 0.0) {
                *err = std::abs(c / sum);
            } else {
                *err = std::abs(c);
            }

            if (*err != *err) {
                /* nan */
                *err = 1.0;
            }

            return (sum);
        }

    } // namespace detail

    XSF_HOST_DEVICE inline double hyperg(double a, double b, double x) {
        double asum, psum, acanc, pcanc, temp;

        /* See if a Kummer transformation will help */
        temp = b - a;
        if (std::abs(temp) < 0.001 * std::abs(a))
            return (exp(x) * hyperg(temp, b, -x));

        /* Try power & asymptotic series, starting from the one that is likely OK */
        if (std::abs(x) < 10 + std::abs(a) + std::abs(b)) {
            psum = detail::hy1f1p(a, b, x, &pcanc);
            if (pcanc < 1.0e-15)
                goto done;
            asum = detail::hy1f1a(a, b, x, &acanc);
        } else {
            psum = detail::hy1f1a(a, b, x, &pcanc);
            if (pcanc < 1.0e-15)
                goto done;
            asum = detail::hy1f1p(a, b, x, &acanc);
        }

        /* Pick the result with less estimated error */

        if (acanc < pcanc) {
            pcanc = acanc;
            psum = asum;
        }

    done:
        if (pcanc > 1.0e-12)
            set_error("hyperg", SF_ERROR_LOSS, NULL);

        return (psum);
    }

} // namespace cephes
} // namespace xsf
