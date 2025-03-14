/* Translated into C++ by SciPy developers in 2024.
 *
 * Original author: Josh Wilson, 2020.
 */

/*
 * Implement sin(pi * x) and cos(pi * x) for real x. Since the periods
 * of these functions are integral (and thus representable in double
 * precision), it's possible to compute them with greater accuracy
 * than sin(x) and cos(x).
 */
#pragma once

#include "../config.h"

namespace xsf {
namespace cephes {

    /* Compute sin(pi * x). */
    template <typename T>
    XSF_HOST_DEVICE T sinpi(T x) {
        T s = 1.0;

        if (x < 0.0) {
            x = -x;
            s = -1.0;
        }

        T r = std::fmod(x, 2.0);
        if (r < 0.5) {
            return s * std::sin(M_PI * r);
        } else if (r > 1.5) {
            return s * std::sin(M_PI * (r - 2.0));
        } else {
            return -s * std::sin(M_PI * (r - 1.0));
        }
    }

    /* Compute cos(pi * x) */
    template <typename T>
    XSF_HOST_DEVICE T cospi(T x) {
        if (x < 0.0) {
            x = -x;
        }

        T r = std::fmod(x, 2.0);
        if (r == 0.5) {
            // We don't want to return -0.0
            return 0.0;
        }
        if (r < 1.0) {
            return -std::sin(M_PI * (r - 0.5));
        } else {
            return std::sin(M_PI * (r - 1.5));
        }
    }
} // namespace cephes
} // namespace xsf
