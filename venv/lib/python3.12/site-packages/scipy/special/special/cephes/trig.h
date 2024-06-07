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

namespace special {
namespace cephes {

    /* Compute sin(pi * x). */
    SPECFUN_HOST_DEVICE double sinpi(double x) {
        double s = 1.0;

        if (x < 0.0) {
            x = -x;
            s = -1.0;
        }

        double r = fmod(x, 2.0);
        if (r < 0.5) {
            return s * sin(M_PI * r);
        } else if (r > 1.5) {
            return s * sin(M_PI * (r - 2.0));
        } else {
            return -s * sin(M_PI * (r - 1.0));
        }
    }

    /* Compute cos(pi * x) */
    SPECFUN_HOST_DEVICE double cospi(double x) {
        if (x < 0.0) {
            x = -x;
        }

        double r = fmod(x, 2.0);
        if (r == 0.5) {
            // We don't want to return -0.0
            return 0.0;
        }
        if (r < 1.0) {
            return -sin(M_PI * (r - 0.5));
        } else {
            return sin(M_PI * (r - 1.5));
        }
    }
} // namespace cephes
} // namespace special
