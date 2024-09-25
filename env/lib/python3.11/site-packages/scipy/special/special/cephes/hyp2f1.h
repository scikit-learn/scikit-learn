/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                      hyp2f1.c
 *
 *      Gauss hypergeometric function   F
 *                                     2 1
 *
 *
 * SYNOPSIS:
 *
 * double a, b, c, x, y, hyp2f1();
 *
 * y = hyp2f1( a, b, c, x );
 *
 *
 * DESCRIPTION:
 *
 *
 *  hyp2f1( a, b, c, x )  =   F ( a, b; c; x )
 *                           2 1
 *
 *           inf.
 *            -   a(a+1)...(a+k) b(b+1)...(b+k)   k+1
 *   =  1 +   >   -----------------------------  x   .
 *            -         c(c+1)...(c+k) (k+1)!
 *          k = 0
 *
 *  Cases addressed are
 *      Tests and escapes for negative integer a, b, or c
 *      Linear transformation if c - a or c - b negative integer
 *      Special case c = a or c = b
 *      Linear transformation for  x near +1
 *      Transformation for x < -0.5
 *      Psi function expansion if x > 0.5 and c - a - b integer
 *      Conditionally, a recurrence on c to make c-a-b > 0
 *
 *      x < -1  AMS 15.3.7 transformation applied (Travis Oliphant)
 *         valid for b,a,c,(b-a) != integer and (c-a),(c-b) != negative integer
 *
 * x >= 1 is rejected (unless special cases are present)
 *
 * The parameters a, b, c are considered to be integer
 * valued if they are within 1.0e-14 of the nearest integer
 * (1.0e-13 for IEEE arithmetic).
 *
 * ACCURACY:
 *
 *
 *               Relative error (-1 < x < 1):
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -1,7        230000      1.2e-11     5.2e-14
 *
 * Several special cases also tested with a, b, c in
 * the range -7 to 7.
 *
 * ERROR MESSAGES:
 *
 * A "partial loss of precision" message is printed if
 * the internally estimated relative error exceeds 1^-12.
 * A "singularity" message is printed on overflow or
 * in cases not addressed (such as x < -1).
 */

/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
 */

#pragma once

#include "../config.h"
#include "../error.h"

#include "const.h"
#include "gamma.h"
#include "psi.h"

namespace special {
namespace cephes {

    namespace detail {
        constexpr double hyp2f1_EPS = 1.0e-13;

        constexpr double hyp2f1_ETHRESH = 1.0e-12;
        constexpr std::uint64_t hyp2f1_MAXITER = 10000;

        /* hys2f1 and hyp2f1ra depend on each other, so we need this prototype */
        SPECFUN_HOST_DEVICE double hyp2f1ra(double a, double b, double c, double x, double *loss);

        /* Defining power series expansion of Gauss hypergeometric function */
        /* The `loss` parameter estimates loss of significance */
        SPECFUN_HOST_DEVICE double hys2f1(double a, double b, double c, double x, double *loss) {
            double f, g, h, k, m, s, u, umax;
            std::uint64_t i;
            int ib, intflag = 0;

            if (std::abs(b) > std::abs(a)) {
                /* Ensure that |a| > |b| ... */
                f = b;
                b = a;
                a = f;
            }

            ib = std::round(b);

            if (std::abs(b - ib) < hyp2f1_EPS && ib <= 0 && std::abs(b) < std::abs(a)) {
                /* .. except when `b` is a smaller negative integer */
                f = b;
                b = a;
                a = f;
                intflag = 1;
            }

            if ((std::abs(a) > std::abs(c) + 1 || intflag) && std::abs(c - a) > 2 && std::abs(a) > 2) {
                /* |a| >> |c| implies that large cancellation error is to be expected.
                 *
                 * We try to reduce it with the recurrence relations
                 */
                return hyp2f1ra(a, b, c, x, loss);
            }

            i = 0;
            umax = 0.0;
            f = a;
            g = b;
            h = c;
            s = 1.0;
            u = 1.0;
            k = 0.0;
            do {
                if (std::abs(h) < hyp2f1_EPS) {
                    *loss = 1.0;
                    return std::numeric_limits<double>::infinity();
                }
                m = k + 1.0;
                u = u * ((f + k) * (g + k) * x / ((h + k) * m));
                s += u;
                k = std::abs(u); /* remember largest term summed */
                if (k > umax)
                    umax = k;
                k = m;
                if (++i > hyp2f1_MAXITER) { /* should never happen */
                    *loss = 1.0;
                    return (s);
                }
            } while (s == 0 || std::abs(u / s) > MACHEP);

            /* return estimated relative error */
            *loss = (MACHEP * umax) / fabs(s) + (MACHEP * i);

            return (s);
        }

        /* Apply transformations for |x| near 1 then call the power series */
        SPECFUN_HOST_DEVICE double hyt2f1(double a, double b, double c, double x, double *loss) {
            double p, q, r, s, t, y, w, d, err, err1;
            double ax, id, d1, d2, e, y1;
            int i, aid, sign;

            int ia, ib, neg_int_a = 0, neg_int_b = 0;

            ia = std::round(a);
            ib = std::round(b);

            if (a <= 0 && std::abs(a - ia) < hyp2f1_EPS) { /* a is a negative integer */
                neg_int_a = 1;
            }

            if (b <= 0 && std::abs(b - ib) < hyp2f1_EPS) { /* b is a negative integer */
                neg_int_b = 1;
            }

            err = 0.0;
            s = 1.0 - x;
            if (x < -0.5 && !(neg_int_a || neg_int_b)) {
                if (b > a)
                    y = std::pow(s, -a) * hys2f1(a, c - b, c, -x / s, &err);

                else
                    y = std::pow(s, -b) * hys2f1(c - a, b, c, -x / s, &err);

                goto done;
            }

            d = c - a - b;
            id = std::round(d); /* nearest integer to d */

            if (x > 0.9 && !(neg_int_a || neg_int_b)) {
                if (std::abs(d - id) > MACHEP) {
                    int sgngam;

                    /* test for integer c-a-b */
                    /* Try the power series first */
                    y = hys2f1(a, b, c, x, &err);
                    if (err < hyp2f1_ETHRESH) {
                        goto done;
                    }
                    /* If power series fails, then apply AMS55 #15.3.6 */
                    q = hys2f1(a, b, 1.0 - d, s, &err);
                    sign = 1;
                    w = lgam_sgn(d, &sgngam);
                    sign *= sgngam;
                    w -= lgam_sgn(c - a, &sgngam);
                    sign *= sgngam;
                    w -= lgam_sgn(c - b, &sgngam);
                    sign *= sgngam;
                    q *= sign * std::exp(w);
                    r = std::pow(s, d) * hys2f1(c - a, c - b, d + 1.0, s, &err1);
                    sign = 1;
                    w = lgam_sgn(-d, &sgngam);
                    sign *= sgngam;
                    w -= lgam_sgn(a, &sgngam);
                    sign *= sgngam;
                    w -= lgam_sgn(b, &sgngam);
                    sign *= sgngam;
                    r *= sign * std::exp(w);
                    y = q + r;

                    q = std::abs(q); /* estimate cancellation error */
                    r = std::abs(r);
                    if (q > r) {
                        r = q;
                    }
                    err += err1 + (MACHEP * r) / y;

                    y *= special::cephes::Gamma(c);
                    goto done;
                } else {
                    /* Psi function expansion, AMS55 #15.3.10, #15.3.11, #15.3.12
                     *
                     * Although AMS55 does not explicitly state it, this expansion fails
                     * for negative integer a or b, since the psi and Gamma functions
                     * involved have poles.
                     */

                    if (id >= 0.0) {
                        e = d;
                        d1 = d;
                        d2 = 0.0;
                        aid = id;
                    } else {
                        e = -d;
                        d1 = 0.0;
                        d2 = d;
                        aid = -id;
                    }

                    ax = std::log(s);

                    /* sum for t = 0 */
                    y = special::cephes::psi(1.0) + special::cephes::psi(1.0 + e) - special::cephes::psi(a + d1) -
                        special::cephes::psi(b + d1) - ax;
                    y /= special::cephes::Gamma(e + 1.0);

                    p = (a + d1) * (b + d1) * s / special::cephes::Gamma(e + 2.0); /* Poch for t=1 */
                    t = 1.0;
                    do {
                        r = special::cephes::psi(1.0 + t) + special::cephes::psi(1.0 + t + e) -
                            special::cephes::psi(a + t + d1) - special::cephes::psi(b + t + d1) - ax;
                        q = p * r;
                        y += q;
                        p *= s * (a + t + d1) / (t + 1.0);
                        p *= (b + t + d1) / (t + 1.0 + e);
                        t += 1.0;
                        if (t > hyp2f1_MAXITER) { /* should never happen */
                            set_error("hyp2f1", SF_ERROR_SLOW, NULL);
                            *loss = 1.0;
                            return std::numeric_limits<double>::quiet_NaN();
                        }
                    } while (y == 0 || std::abs(q / y) > hyp2f1_EPS);

                    if (id == 0.0) {
                        y *= special::cephes::Gamma(c) / (special::cephes::Gamma(a) * special::cephes::Gamma(b));
                        goto psidon;
                    }

                    y1 = 1.0;

                    if (aid == 1)
                        goto nosum;

                    t = 0.0;
                    p = 1.0;
                    for (i = 1; i < aid; i++) {
                        r = 1.0 - e + t;
                        p *= s * (a + t + d2) * (b + t + d2) / r;
                        t += 1.0;
                        p /= t;
                        y1 += p;
                    }
                nosum:
                    p = special::cephes::Gamma(c);
                    y1 *= special::cephes::Gamma(e) * p /
                          (special::cephes::Gamma(a + d1) * special::cephes::Gamma(b + d1));

                    y *= p / (special::cephes::Gamma(a + d2) * special::cephes::Gamma(b + d2));
                    if ((aid & 1) != 0)
                        y = -y;

                    q = std::pow(s, id); /* s to the id power */
                    if (id > 0.0)
                        y *= q;
                    else
                        y1 *= q;

                    y += y1;
                psidon:
                    goto done;
                }
            }

            /* Use defining power series if no special cases */
            y = hys2f1(a, b, c, x, &err);

        done:
            *loss = err;
            return (y);
        }

        /*
          15.4.2 Abramowitz & Stegun.
        */
        SPECFUN_HOST_DEVICE double hyp2f1_neg_c_equal_bc(double a, double b, double x) {
            double k;
            double collector = 1;
            double sum = 1;
            double collector_max = 1;

            if (!(std::abs(b) < 1e5)) {
                return std::numeric_limits<double>::quiet_NaN();
            }

            for (k = 1; k <= -b; k++) {
                collector *= (a + k - 1) * x / k;
                collector_max = std::fmax(std::abs(collector), collector_max);
                sum += collector;
            }

            if (1e-16 * (1 + collector_max / std::abs(sum)) > 1e-7) {
                return std::numeric_limits<double>::quiet_NaN();
            }

            return sum;
        }

        /*
         * Evaluate hypergeometric function by two-term recurrence in `a`.
         *
         * This avoids some of the loss of precision in the strongly alternating
         * hypergeometric series, and can be used to reduce the `a` and `b` parameters
         * to smaller values.
         *
         * AMS55 #15.2.10
         */
        SPECFUN_HOST_DEVICE double hyp2f1ra(double a, double b, double c, double x, double *loss) {
            double f2, f1, f0;
            int n;
            double t, err, da;

            /* Don't cross c or zero */
            if ((c < 0 && a <= c) || (c >= 0 && a >= c)) {
                da = std::round(a - c);
            } else {
                da = std::round(a);
            }
            t = a - da;

            *loss = 0;

            SPECFUN_ASSERT(da != 0);

            if (std::abs(da) > hyp2f1_MAXITER) {
                /* Too expensive to compute this value, so give up */
                set_error("hyp2f1", SF_ERROR_NO_RESULT, NULL);
                *loss = 1.0;
                return std::numeric_limits<double>::quiet_NaN();
            }

            if (da < 0) {
                /* Recurse down */
                f2 = 0;
                f1 = hys2f1(t, b, c, x, &err);
                *loss += err;
                f0 = hys2f1(t - 1, b, c, x, &err);
                *loss += err;
                t -= 1;
                for (n = 1; n < -da; ++n) {
                    f2 = f1;
                    f1 = f0;
                    f0 = -(2 * t - c - t * x + b * x) / (c - t) * f1 - t * (x - 1) / (c - t) * f2;
                    t -= 1;
                }
            } else {
                /* Recurse up */
                f2 = 0;
                f1 = hys2f1(t, b, c, x, &err);
                *loss += err;
                f0 = hys2f1(t + 1, b, c, x, &err);
                *loss += err;
                t += 1;
                for (n = 1; n < da; ++n) {
                    f2 = f1;
                    f1 = f0;
                    f0 = -((2 * t - c - t * x + b * x) * f1 + (c - t) * f2) / (t * (x - 1));
                    t += 1;
                }
            }

            return f0;
        }
    } // namespace detail

    SPECFUN_HOST_DEVICE double hyp2f1(double a, double b, double c, double x) {
        double d, d1, d2, e;
        double p, q, r, s, y, ax;
        double ia, ib, ic, id, err;
        double t1;
        int i, aid;
        int neg_int_a = 0, neg_int_b = 0;
        int neg_int_ca_or_cb = 0;

        err = 0.0;
        ax = std::abs(x);
        s = 1.0 - x;
        ia = std::round(a); /* nearest integer to a */
        ib = std::round(b);

        if (x == 0.0) {
            return 1.0;
        }

        d = c - a - b;
        id = std::round(d);

        if ((a == 0 || b == 0) && c != 0) {
            return 1.0;
        }

        if (a <= 0 && std::abs(a - ia) < detail::hyp2f1_EPS) { /* a is a negative integer */
            neg_int_a = 1;
        }

        if (b <= 0 && std::abs(b - ib) < detail::hyp2f1_EPS) { /* b is a negative integer */
            neg_int_b = 1;
        }

        if (d <= -1 && !(std::abs(d - id) > detail::hyp2f1_EPS && s < 0) && !(neg_int_a || neg_int_b)) {
            return std::pow(s, d) * hyp2f1(c - a, c - b, c, x);
        }
        if (d <= 0 && x == 1 && !(neg_int_a || neg_int_b))
            goto hypdiv;

        if (ax < 1.0 || x == -1.0) {
            /* 2F1(a,b;b;x) = (1-x)**(-a) */
            if (std::abs(b - c) < detail::hyp2f1_EPS) { /* b = c */
                if (neg_int_b) {
                    y = detail::hyp2f1_neg_c_equal_bc(a, b, x);
                } else {
                    y = std::pow(s, -a); /* s to the -a power */
                }
                goto hypdon;
            }
            if (std::abs(a - c) < detail::hyp2f1_EPS) { /* a = c */
                y = std::pow(s, -b);                    /* s to the -b power */
                goto hypdon;
            }
        }

        if (c <= 0.0) {
            ic = std::round(c);                          /* nearest integer to c */
            if (std::abs(c - ic) < detail::hyp2f1_EPS) { /* c is a negative integer */
                /* check if termination before explosion */
                if (neg_int_a && (ia > ic))
                    goto hypok;
                if (neg_int_b && (ib > ic))
                    goto hypok;
                goto hypdiv;
            }
        }

        if (neg_int_a || neg_int_b) /* function is a polynomial */
            goto hypok;

        t1 = std::abs(b - a);
        if (x < -2.0 && std::abs(t1 - round(t1)) > detail::hyp2f1_EPS) {
            /* This transform has a pole for b-a integer, and
             * may produce large cancellation errors for |1/x| close 1
             */
            p = hyp2f1(a, 1 - c + a, 1 - b + a, 1.0 / x);
            q = hyp2f1(b, 1 - c + b, 1 - a + b, 1.0 / x);
            p *= std::pow(-x, -a);
            q *= std::pow(-x, -b);
            t1 = Gamma(c);
            s = t1 * Gamma(b - a) / (Gamma(b) * Gamma(c - a));
            y = t1 * Gamma(a - b) / (Gamma(a) * Gamma(c - b));
            return s * p + y * q;
        } else if (x < -1.0) {
            if (std::abs(a) < std::abs(b)) {
                return std::pow(s, -a) * hyp2f1(a, c - b, c, x / (x - 1));
            } else {
                return std::pow(s, -b) * hyp2f1(b, c - a, c, x / (x - 1));
            }
        }

        if (ax > 1.0) /* series diverges  */
            goto hypdiv;

        p = c - a;
        ia = std::round(p);                                         /* nearest integer to c-a */
        if ((ia <= 0.0) && (std::abs(p - ia) < detail::hyp2f1_EPS)) /* negative int c - a */
            neg_int_ca_or_cb = 1;

        r = c - b;
        ib = std::round(r);                                         /* nearest integer to c-b */
        if ((ib <= 0.0) && (std::abs(r - ib) < detail::hyp2f1_EPS)) /* negative int c - b */
            neg_int_ca_or_cb = 1;

        id = std::round(d); /* nearest integer to d */
        q = std::abs(d - id);

        /* Thanks to Christian Burger <BURGER@DMRHRZ11.HRZ.Uni-Marburg.DE>
         * for reporting a bug here.  */
        if (std::abs(ax - 1.0) < detail::hyp2f1_EPS) { /* |x| == 1.0   */
            if (x > 0.0) {
                if (neg_int_ca_or_cb) {
                    if (d >= 0.0)
                        goto hypf;
                    else
                        goto hypdiv;
                }
                if (d <= 0.0)
                    goto hypdiv;
                y = Gamma(c) * Gamma(d) / (Gamma(p) * Gamma(r));
                goto hypdon;
            }
            if (d <= -1.0)
                goto hypdiv;
        }

        /* Conditionally make d > 0 by recurrence on c
         * AMS55 #15.2.27
         */
        if (d < 0.0) {
            /* Try the power series first */
            y = detail::hyt2f1(a, b, c, x, &err);
            if (err < detail::hyp2f1_ETHRESH)
                goto hypdon;
            /* Apply the recurrence if power series fails */
            err = 0.0;
            aid = 2 - id;
            e = c + aid;
            d2 = hyp2f1(a, b, e, x);
            d1 = hyp2f1(a, b, e + 1.0, x);
            q = a + b + 1.0;
            for (i = 0; i < aid; i++) {
                r = e - 1.0;
                y = (e * (r - (2.0 * e - q) * x) * d2 + (e - a) * (e - b) * x * d1) / (e * r * s);
                e = r;
                d1 = d2;
                d2 = y;
            }
            goto hypdon;
        }

        if (neg_int_ca_or_cb) {
            goto hypf; /* negative integer c-a or c-b */
        }

    hypok:
        y = detail::hyt2f1(a, b, c, x, &err);

    hypdon:
        if (err > detail::hyp2f1_ETHRESH) {
            set_error("hyp2f1", SF_ERROR_LOSS, NULL);
            /*      printf( "Estimated err = %.2e\n", err ); */
        }
        return (y);

        /* The transformation for c-a or c-b negative integer
         * AMS55 #15.3.3
         */
    hypf:
        y = std::pow(s, d) * detail::hys2f1(c - a, c - b, c, x, &err);
        goto hypdon;

        /* The alarm exit */
    hypdiv:
        set_error("hyp2f1", SF_ERROR_OVERFLOW, NULL);
        return std::numeric_limits<double>::infinity();
    }

} // namespace cephes
} // namespace special
