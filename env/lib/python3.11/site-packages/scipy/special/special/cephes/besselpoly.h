/* Translated into C++ by SciPy developers in 2024.
 *
 * This was not part of the original cephes library.
 */
#pragma once

#include "../config.h"
#include "gamma.h"

namespace special {
namespace cephes {
    namespace detail {

        constexpr double besselpoly_EPS = 1.0e-17;
    }

    SPECFUN_HOST_DEVICE inline double besselpoly(double a, double lambda, double nu) {

        int m, factor = 0;
        double Sm, relerr, Sol;
        double sum = 0.0;

        /* Special handling for a = 0.0 */
        if (a == 0.0) {
            if (nu == 0.0) {
                return 1.0 / (lambda + 1);
            } else {
                return 0.0;
            }
        }
        /* Special handling for negative and integer nu */
        if ((nu < 0) && (std::floor(nu) == nu)) {
            nu = -nu;
            factor = static_cast<int>(nu) % 2;
        }
        Sm = std::exp(nu * std::log(a)) / (Gamma(nu + 1) * (lambda + nu + 1));
        m = 0;
        do {
            sum += Sm;
            Sol = Sm;
            Sm *= -a * a * (lambda + nu + 1 + 2 * m) / ((nu + m + 1) * (m + 1) * (lambda + nu + 1 + 2 * m + 2));
            m++;
            relerr = std::abs((Sm - Sol) / Sm);
        } while (relerr > detail::besselpoly_EPS && m < 1000);
        if (!factor)
            return sum;
        else
            return -sum;
    }
} // namespace cephes
} // namespace special
