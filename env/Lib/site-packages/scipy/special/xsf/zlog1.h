/* Translated from Cython into C++ by SciPy developers in 2023.
 *
 * Original author: Josh Wilson, 2016.
 */

#pragma once

#include "config.h"

namespace xsf {
namespace detail {

    XSF_HOST_DEVICE inline std::complex<double> zlog1(std::complex<double> z) {
        /* Compute log, paying special attention to accuracy around 1. We
         * implement this ourselves because some systems (most notably the
         * Travis CI machines) are weak in this regime. */
        std::complex<double> coeff = -1.0;
        std::complex<double> res = 0.0;

        if (std::abs(z - 1.0) > 0.1) {
            return std::log(z);
        }

        z -= 1.0;
        for (int n = 1; n < 17; n++) {
            coeff *= -z;
            res += coeff / static_cast<double>(n);
            if (std::abs(res / coeff) < std::numeric_limits<double>::epsilon()) {
                break;
            }
        }
        return res;
    }
} // namespace detail
} // namespace xsf
