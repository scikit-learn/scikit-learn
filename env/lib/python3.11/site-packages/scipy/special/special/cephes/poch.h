/*
 * Pochhammer symbol (a)_m = gamma(a + m) / gamma(a)
 */

#pragma once

#include "../config.h"
#include "gamma.h"

namespace special {
namespace cephes {

    namespace detail {

        SPECFUN_HOST_DEVICE inline double is_nonpos_int(double x) {
            return x <= 0 && x == std::ceil(x) && std::abs(x) < 1e13;
        }
    } // namespace detail

    SPECFUN_HOST_DEVICE inline double poch(double a, double m) {
        double r = 1.0;

        /*
         * 1. Reduce magnitude of `m` to |m| < 1 by using recurrence relations.
         *
         * This may end up in over/underflow, but then the function itself either
         * diverges or goes to zero. In case the remainder goes to the opposite
         * direction, we end up returning 0*INF = NAN, which is OK.
         */

        /* Recurse down */
        while (m >= 1.0) {
            if (a + m == 1) {
                break;
            }
            m -= 1.0;
            r *= (a + m);
            if (!std::isfinite(r) || r == 0) {
                break;
            }
        }

        /* Recurse up */
        while (m <= -1.0) {
            if (a + m == 0) {
                break;
            }
            r /= (a + m);
            m += 1.0;
            if (!std::isfinite(r) || r == 0) {
                break;
            }
        }

        /*
         * 2. Evaluate function with reduced `m`
         *
         * Now either `m` is not big, or the `r` product has over/underflown.
         * If so, the function itself does similarly.
         */

        if (m == 0) {
            /* Easy case */
            return r;
        } else if (a > 1e4 && std::abs(m) <= 1) {
            /* Avoid loss of precision */
            return r * std::pow(a, m) *
                   (1 + m * (m - 1) / (2 * a) + m * (m - 1) * (m - 2) * (3 * m - 1) / (24 * a * a) +
                    m * m * (m - 1) * (m - 1) * (m - 2) * (m - 3) / (48 * a * a * a));
        }

        /* Check for infinity */
        if (detail::is_nonpos_int(a + m) && !detail::is_nonpos_int(a) && a + m != m) {
            return std::numeric_limits<double>::infinity();
        }

        /* Check for zero */
        if (!detail::is_nonpos_int(a + m) && detail::is_nonpos_int(a)) {
            return 0;
        }

        return r * std::exp(lgam(a + m) - lgam(a)) * gammasgn(a + m) * gammasgn(a);
    }
} // namespace cephes
} // namespace special
