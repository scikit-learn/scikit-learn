/* Translated into C++ by SciPy developers in 2024. */

/*                                                     unity.c
 *
 * Relative error approximations for function arguments near
 * unity.
 *
 *    log1p(x) = log(1+x)
 *    expm1(x) = exp(x) - 1
 *    cosm1(x) = cos(x) - 1
 *    lgam1p(x) = lgam(1+x)
 *
 */

/* Scipy changes:
 * - 06-10-2016: added lgam1p
 */
#pragma once

#include "../config.h"

#include "const.h"
#include "gamma.h"
#include "polevl.h"
#include "zeta.h"

namespace xsf {
namespace cephes {

    namespace detail {

        /* log1p(x) = log(1 + x)  */

        /* Coefficients for log(1+x) = x - x**2/2 + x**3 P(x)/Q(x)
         * 1/sqrt(2) <= x < sqrt(2)
         * Theoretical peak relative error = 2.32e-20
         */

        constexpr double unity_LP[] = {
            4.5270000862445199635215E-5, 4.9854102823193375972212E-1, 6.5787325942061044846969E0,
            2.9911919328553073277375E1,  6.0949667980987787057556E1,  5.7112963590585538103336E1,
            2.0039553499201281259648E1,
        };

        constexpr double unity_LQ[] = {
            /* 1.0000000000000000000000E0, */
            1.5062909083469192043167E1, 8.3047565967967209469434E1, 2.2176239823732856465394E2,
            3.0909872225312059774938E2, 2.1642788614495947685003E2, 6.0118660497603843919306E1,
        };

    } // namespace detail

    XSF_HOST_DEVICE inline double log1p(double x) {
        double z;

        z = 1.0 + x;
        if ((z < M_SQRT1_2) || (z > M_SQRT2))
            return (std::log(z));
        z = x * x;
        z = -0.5 * z + x * (z * polevl(x, detail::unity_LP, 6) / p1evl(x, detail::unity_LQ, 6));
        return (x + z);
    }

    /* log(1 + x) - x */
    XSF_HOST_DEVICE inline double log1pmx(double x) {
        if (std::abs(x) < 0.5) {
            uint64_t n;
            double xfac = x;
            double term;
            double res = 0;

            for (n = 2; n < detail::MAXITER; n++) {
                xfac *= -x;
                term = xfac / n;
                res += term;
                if (std::abs(term) < detail::MACHEP * std::abs(res)) {
                    break;
                }
            }
            return res;
        } else {
            return log1p(x) - x;
        }
    }

    /* expm1(x) = exp(x) - 1  */

    /*  e^x =  1 + 2x P(x^2)/( Q(x^2) - P(x^2) )
     * -0.5 <= x <= 0.5
     */

    namespace detail {

        constexpr double unity_EP[3] = {
            1.2617719307481059087798E-4,
            3.0299440770744196129956E-2,
            9.9999999999999999991025E-1,
        };

        constexpr double unity_EQ[4] = {
            3.0019850513866445504159E-6,
            2.5244834034968410419224E-3,
            2.2726554820815502876593E-1,
            2.0000000000000000000897E0,
        };

    } // namespace detail

    XSF_HOST_DEVICE inline double expm1(double x) {
        double r, xx;

        if (!std::isfinite(x)) {
            if (std::isnan(x)) {
                return x;
            } else if (x > 0) {
                return x;
            } else {
                return -1.0;
            }
        }
        if ((x < -0.5) || (x > 0.5))
            return (std::exp(x) - 1.0);
        xx = x * x;
        r = x * polevl(xx, detail::unity_EP, 2);
        r = r / (polevl(xx, detail::unity_EQ, 3) - r);
        return (r + r);
    }

    /* cosm1(x) = cos(x) - 1  */

    namespace detail {
        constexpr double unity_coscof[7] = {
            4.7377507964246204691685E-14, -1.1470284843425359765671E-11, 2.0876754287081521758361E-9,
            -2.7557319214999787979814E-7, 2.4801587301570552304991E-5,   -1.3888888888888872993737E-3,
            4.1666666666666666609054E-2,
        };

    }

    XSF_HOST_DEVICE inline double cosm1(double x) {
        double xx;

        if ((x < -M_PI_4) || (x > M_PI_4))
            return (std::cos(x) - 1.0);
        xx = x * x;
        xx = -0.5 * xx + xx * xx * polevl(xx, detail::unity_coscof, 6);
        return xx;
    }

    namespace detail {
        /* Compute lgam(x + 1) around x = 0 using its Taylor series. */
        XSF_HOST_DEVICE inline double lgam1p_taylor(double x) {
            int n;
            double xfac, coeff, res;

            if (x == 0) {
                return 0;
            }
            res = -SCIPY_EULER * x;
            xfac = -x;
            for (n = 2; n < 42; n++) {
                xfac *= -x;
                coeff = xsf::cephes::zeta(n, 1) * xfac / n;
                res += coeff;
                if (std::abs(coeff) < detail::MACHEP * std::abs(res)) {
                    break;
                }
            }

            return res;
        }
    } // namespace detail

    /* Compute lgam(x + 1). */
    XSF_HOST_DEVICE inline double lgam1p(double x) {
        if (std::abs(x) <= 0.5) {
            return detail::lgam1p_taylor(x);
        } else if (std::abs(x - 1) < 0.5) {
            return std::log(x) + detail::lgam1p_taylor(x - 1);
        } else {
            return lgam(x + 1);
        }
    }

} // namespace cephes
} // namespace xsf
