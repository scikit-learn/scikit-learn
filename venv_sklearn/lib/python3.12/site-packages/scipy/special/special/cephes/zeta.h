/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     zeta.c
 *
 *     Riemann zeta function of two arguments
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, q, y, zeta();
 *
 * y = zeta( x, q );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *
 *                 inf.
 *                  -        -x
 *   zeta(x,q)  =   >   (k+q)
 *                  -
 *                 k=0
 *
 * where x > 1 and q is not a negative integer or zero.
 * The Euler-Maclaurin summation formula is used to obtain
 * the expansion
 *
 *                n
 *                -       -x
 * zeta(x,q)  =   >  (k+q)
 *                -
 *               k=1
 *
 *           1-x                 inf.  B   x(x+1)...(x+2j)
 *      (n+q)           1         -     2j
 *  +  ---------  -  -------  +   >    --------------------
 *        x-1              x      -                   x+2j+1
 *                   2(n+q)      j=1       (2j)! (n+q)
 *
 * where the B2j are Bernoulli numbers.  Note that (see zetac.c)
 * zeta(x,1) = zetac(x) + 1.
 *
 *
 *
 * ACCURACY:
 *
 *
 *
 * REFERENCE:
 *
 * Gradshteyn, I. S., and I. M. Ryzhik, Tables of Integrals,
 * Series, and Products, p. 1073; Academic Press, 1980.
 *
 */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
#pragma once

#include "../config.h"
#include "../error.h"
#include "const.h"

namespace special {
namespace cephes {

    namespace detail {
        /* Expansion coefficients
         * for Euler-Maclaurin summation formula
         * (2k)! / B2k
         * where B2k are Bernoulli numbers
         */
        constexpr double zeta_A[] = {
            12.0,
            -720.0,
            30240.0,
            -1209600.0,
            47900160.0,
            -1.8924375803183791606e9, /*1.307674368e12/691 */
            7.47242496e10,
            -2.950130727918164224e12,  /*1.067062284288e16/3617 */
            1.1646782814350067249e14,  /*5.109094217170944e18/43867 */
            -4.5979787224074726105e15, /*8.028576626982912e20/174611 */
            1.8152105401943546773e17,  /*1.5511210043330985984e23/854513 */
            -7.1661652561756670113e18  /*1.6938241367317436694528e27/236364091 */
        };

        /* 30 Nov 86 -- error in third coefficient fixed */
    } // namespace detail

    SPECFUN_HOST_DEVICE double inline zeta(double x, double q) {
        int i;
        double a, b, k, s, t, w;

        if (x == 1.0)
            goto retinf;

        if (x < 1.0) {
        domerr:
            set_error("zeta", SF_ERROR_DOMAIN, NULL);
            return (std::numeric_limits<double>::quiet_NaN());
        }

        if (q <= 0.0) {
            if (q == floor(q)) {
                set_error("zeta", SF_ERROR_SINGULAR, NULL);
            retinf:
                return (std::numeric_limits<double>::infinity());
            }
            if (x != std::floor(x))
                goto domerr; /* because q^-x not defined */
        }

        /* Asymptotic expansion
         * https://dlmf.nist.gov/25.11#E43
         */
        if (q > 1e8) {
            return (1 / (x - 1) + 1 / (2 * q)) * std::pow(q, 1 - x);
        }

        /* Euler-Maclaurin summation formula */

        /* Permit negative q but continue sum until n+q > +9 .
         * This case should be handled by a reflection formula.
         * If q<0 and x is an integer, there is a relation to
         * the polyGamma function.
         */
        s = std::pow(q, -x);
        a = q;
        i = 0;
        b = 0.0;
        while ((i < 9) || (a <= 9.0)) {
            i += 1;
            a += 1.0;
            b = std::pow(a, -x);
            s += b;
            if (std::abs(b / s) < detail::MACHEP)
                goto done;
        }

        w = a;
        s += b * w / (x - 1.0);
        s -= 0.5 * b;
        a = 1.0;
        k = 0.0;
        for (i = 0; i < 12; i++) {
            a *= x + k;
            b /= w;
            t = a * b / detail::zeta_A[i];
            s = s + t;
            t = std::abs(t / s);
            if (t < detail::MACHEP)
                goto done;
            k += 1.0;
            a *= x + k;
            b /= w;
            k += 1.0;
        }
    done:
        return (s);
    }

} // namespace cephes
} // namespace special
