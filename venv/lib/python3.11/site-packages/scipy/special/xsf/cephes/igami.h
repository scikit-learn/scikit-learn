/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*
 * (C) Copyright John Maddock 2006.
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt)
 */
#pragma once

#include "../config.h"
#include "../error.h"

#include "const.h"
#include "gamma.h"
#include "igam.h"
#include "polevl.h"

namespace xsf {
namespace cephes {

    namespace detail {

        XSF_HOST_DEVICE double find_inverse_s(double p, double q) {
            /*
             * Computation of the Incomplete Gamma Function Ratios and their Inverse
             * ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
             * ACM Transactions on Mathematical Software, Vol. 12, No. 4,
             * December 1986, Pages 377-393.
             *
             * See equation 32.
             */
            double s, t;
            constexpr double a[4] = {0.213623493715853, 4.28342155967104, 11.6616720288968, 3.31125922108741};
            constexpr double b[5] = {0.3611708101884203e-1, 1.27364489782223, 6.40691597760039, 6.61053765625462, 1};

            if (p < 0.5) {
                t = std::sqrt(-2 * std::log(p));
            } else {
                t = std::sqrt(-2 * std::log(q));
            }
            s = t - polevl(t, a, 3) / polevl(t, b, 4);
            if (p < 0.5)
                s = -s;
            return s;
        }

        XSF_HOST_DEVICE inline double didonato_SN(double a, double x, unsigned N, double tolerance) {
            /*
             * Computation of the Incomplete Gamma Function Ratios and their Inverse
             * ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
             * ACM Transactions on Mathematical Software, Vol. 12, No. 4,
             * December 1986, Pages 377-393.
             *
             * See equation 34.
             */
            double sum = 1.0;

            if (N >= 1) {
                unsigned i;
                double partial = x / (a + 1);

                sum += partial;
                for (i = 2; i <= N; ++i) {
                    partial *= x / (a + i);
                    sum += partial;
                    if (partial < tolerance) {
                        break;
                    }
                }
            }
            return sum;
        }

        XSF_HOST_DEVICE inline double find_inverse_gamma(double a, double p, double q) {
            /*
             * In order to understand what's going on here, you will
             * need to refer to:
             *
             * Computation of the Incomplete Gamma Function Ratios and their Inverse
             * ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
             * ACM Transactions on Mathematical Software, Vol. 12, No. 4,
             * December 1986, Pages 377-393.
             */
            double result;

            if (a == 1) {
                if (q > 0.9) {
                    result = -std::log1p(-p);
                } else {
                    result = -std::log(q);
                }
            } else if (a < 1) {
                double g = xsf::cephes::Gamma(a);
                double b = q * g;

                if ((b > 0.6) || ((b >= 0.45) && (a >= 0.3))) {
                    /* DiDonato & Morris Eq 21:
                     *
                     * There is a slight variation from DiDonato and Morris here:
                     * the first form given here is unstable when p is close to 1,
                     * making it impossible to compute the inverse of Q(a,x) for small
                     * q. Fortunately the second form works perfectly well in this case.
                     */
                    double u;
                    if ((b * q > 1e-8) && (q > 1e-5)) {
                        u = std::pow(p * g * a, 1 / a);
                    } else {
                        u = std::exp((-q / a) - SCIPY_EULER);
                    }
                    result = u / (1 - (u / (a + 1)));
                } else if ((a < 0.3) && (b >= 0.35)) {
                    /* DiDonato & Morris Eq 22: */
                    double t = std::exp(-SCIPY_EULER - b);
                    double u = t * std::exp(t);
                    result = t * std::exp(u);
                } else if ((b > 0.15) || (a >= 0.3)) {
                    /* DiDonato & Morris Eq 23: */
                    double y = -std::log(b);
                    double u = y - (1 - a) * std::log(y);
                    result = y - (1 - a) * std::log(u) - std::log(1 + (1 - a) / (1 + u));
                } else if (b > 0.1) {
                    /* DiDonato & Morris Eq 24: */
                    double y = -std::log(b);
                    double u = y - (1 - a) * std::log(y);
                    result = y - (1 - a) * std::log(u) -
                             std::log((u * u + 2 * (3 - a) * u + (2 - a) * (3 - a)) / (u * u + (5 - a) * u + 2));
                } else {
                    /* DiDonato & Morris Eq 25: */
                    double y = -std::log(b);
                    double c1 = (a - 1) * std::log(y);
                    double c1_2 = c1 * c1;
                    double c1_3 = c1_2 * c1;
                    double c1_4 = c1_2 * c1_2;
                    double a_2 = a * a;
                    double a_3 = a_2 * a;

                    double c2 = (a - 1) * (1 + c1);
                    double c3 = (a - 1) * (-(c1_2 / 2) + (a - 2) * c1 + (3 * a - 5) / 2);
                    double c4 = (a - 1) * ((c1_3 / 3) - (3 * a - 5) * c1_2 / 2 + (a_2 - 6 * a + 7) * c1 +
                                           (11 * a_2 - 46 * a + 47) / 6);
                    double c5 = (a - 1) * (-(c1_4 / 4) + (11 * a - 17) * c1_3 / 6 + (-3 * a_2 + 13 * a - 13) * c1_2 +
                                           (2 * a_3 - 25 * a_2 + 72 * a - 61) * c1 / 2 +
                                           (25 * a_3 - 195 * a_2 + 477 * a - 379) / 12);

                    double y_2 = y * y;
                    double y_3 = y_2 * y;
                    double y_4 = y_2 * y_2;
                    result = y + c1 + (c2 / y) + (c3 / y_2) + (c4 / y_3) + (c5 / y_4);
                }
            } else {
                /* DiDonato and Morris Eq 31: */
                double s = find_inverse_s(p, q);

                double s_2 = s * s;
                double s_3 = s_2 * s;
                double s_4 = s_2 * s_2;
                double s_5 = s_4 * s;
                double ra = std::sqrt(a);

                double w = a + s * ra + (s_2 - 1) / 3;
                w += (s_3 - 7 * s) / (36 * ra);
                w -= (3 * s_4 + 7 * s_2 - 16) / (810 * a);
                w += (9 * s_5 + 256 * s_3 - 433 * s) / (38880 * a * ra);

                if ((a >= 500) && (std::abs(1 - w / a) < 1e-6)) {
                    result = w;
                } else if (p > 0.5) {
                    if (w < 3 * a) {
                        result = w;
                    } else {
                        double D = std::fmax(2, a * (a - 1));
                        double lg = xsf::cephes::lgam(a);
                        double lb = std::log(q) + lg;
                        if (lb < -D * 2.3) {
                            /* DiDonato and Morris Eq 25: */
                            double y = -lb;
                            double c1 = (a - 1) * std::log(y);
                            double c1_2 = c1 * c1;
                            double c1_3 = c1_2 * c1;
                            double c1_4 = c1_2 * c1_2;
                            double a_2 = a * a;
                            double a_3 = a_2 * a;

                            double c2 = (a - 1) * (1 + c1);
                            double c3 = (a - 1) * (-(c1_2 / 2) + (a - 2) * c1 + (3 * a - 5) / 2);
                            double c4 = (a - 1) * ((c1_3 / 3) - (3 * a - 5) * c1_2 / 2 + (a_2 - 6 * a + 7) * c1 +
                                                   (11 * a_2 - 46 * a + 47) / 6);
                            double c5 =
                                (a - 1) * (-(c1_4 / 4) + (11 * a - 17) * c1_3 / 6 + (-3 * a_2 + 13 * a - 13) * c1_2 +
                                           (2 * a_3 - 25 * a_2 + 72 * a - 61) * c1 / 2 +
                                           (25 * a_3 - 195 * a_2 + 477 * a - 379) / 12);

                            double y_2 = y * y;
                            double y_3 = y_2 * y;
                            double y_4 = y_2 * y_2;
                            result = y + c1 + (c2 / y) + (c3 / y_2) + (c4 / y_3) + (c5 / y_4);
                        } else {
                            /* DiDonato and Morris Eq 33: */
                            double u = -lb + (a - 1) * std::log(w) - std::log(1 + (1 - a) / (1 + w));
                            result = -lb + (a - 1) * std::log(u) - std::log(1 + (1 - a) / (1 + u));
                        }
                    }
                } else {
                    double z = w;
                    double ap1 = a + 1;
                    double ap2 = a + 2;
                    if (w < 0.15 * ap1) {
                        /* DiDonato and Morris Eq 35: */
                        double v = std::log(p) + xsf::cephes::lgam(ap1);
                        z = std::exp((v + w) / a);
                        s = std::log1p(z / ap1 * (1 + z / ap2));
                        z = std::exp((v + z - s) / a);
                        s = std::log1p(z / ap1 * (1 + z / ap2));
                        z = std::exp((v + z - s) / a);
                        s = std::log1p(z / ap1 * (1 + z / ap2 * (1 + z / (a + 3))));
                        z = std::exp((v + z - s) / a);
                    }

                    if ((z <= 0.01 * ap1) || (z > 0.7 * ap1)) {
                        result = z;
                    } else {
                        /* DiDonato and Morris Eq 36: */
                        double ls = std::log(didonato_SN(a, z, 100, 1e-4));
                        double v = std::log(p) + xsf::cephes::lgam(ap1);
                        z = std::exp((v + z - ls) / a);
                        result = z * (1 - (a * std::log(z) - z - v + ls) / (a - z));
                    }
                }
            }
            return result;
        }

    } // namespace detail

    XSF_HOST_DEVICE inline double igamci(double a, double q);

    XSF_HOST_DEVICE inline double igami(double a, double p) {
        int i;
        double x, fac, f_fp, fpp_fp;

        if (std::isnan(a) || std::isnan(p)) {
            return std::numeric_limits<double>::quiet_NaN();
            ;
        } else if ((a < 0) || (p < 0) || (p > 1)) {
            set_error("gammaincinv", SF_ERROR_DOMAIN, NULL);
        } else if (p == 0.0) {
            return 0.0;
        } else if (p == 1.0) {
            return std::numeric_limits<double>::infinity();
        } else if (p > 0.9) {
            return igamci(a, 1 - p);
        }

        x = detail::find_inverse_gamma(a, p, 1 - p);
        /* Halley's method */
        for (i = 0; i < 3; i++) {
            fac = detail::igam_fac(a, x);
            if (fac == 0.0) {
                return x;
            }
            f_fp = (igam(a, x) - p) * x / fac;
            /* The ratio of the first and second derivatives simplifies */
            fpp_fp = -1.0 + (a - 1) / x;
            if (std::isinf(fpp_fp)) {
                /* Resort to Newton's method in the case of overflow */
                x = x - f_fp;
            } else {
                x = x - f_fp / (1.0 - 0.5 * f_fp * fpp_fp);
            }
        }

        return x;
    }

    XSF_HOST_DEVICE inline double igamci(double a, double q) {
        int i;
        double x, fac, f_fp, fpp_fp;

        if (std::isnan(a) || std::isnan(q)) {
            return std::numeric_limits<double>::quiet_NaN();
        } else if ((a < 0.0) || (q < 0.0) || (q > 1.0)) {
            set_error("gammainccinv", SF_ERROR_DOMAIN, NULL);
        } else if (q == 0.0) {
            return std::numeric_limits<double>::infinity();
        } else if (q == 1.0) {
            return 0.0;
        } else if (q > 0.9) {
            return igami(a, 1 - q);
        }

        x = detail::find_inverse_gamma(a, 1 - q, q);
        for (i = 0; i < 3; i++) {
            fac = detail::igam_fac(a, x);
            if (fac == 0.0) {
                return x;
            }
            f_fp = (igamc(a, x) - q) * x / (-fac);
            fpp_fp = -1.0 + (a - 1) / x;
            if (std::isinf(fpp_fp)) {
                x = x - f_fp;
            } else {
                x = x - f_fp / (1.0 - 0.5 * f_fp * fpp_fp);
            }
        }

        return x;
    }

} // namespace cephes
} // namespace xsf
