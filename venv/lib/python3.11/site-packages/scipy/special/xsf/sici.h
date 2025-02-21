/* Translated from Cython into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/* Implementation of sin/cos/sinh/cosh integrals for complex arguments
 *
 * Sources
 * [1] Fredrik Johansson and others. mpmath: a Python library for
 *     arbitrary-precision floating-point arithmetic (version 0.19),
 *     December 2013. http://mpmath.org/.
 * [2] NIST, "Digital Library of Mathematical Functions",
 *     https://dlmf.nist.gov/
 */

#pragma once

#include "config.h"
#include "error.h"

#include "expint.h"
#include "cephes/const.h"
#include "cephes/sici.h"
#include "cephes/shichi.h"

namespace xsf {
namespace detail {
    
    XSF_HOST_DEVICE inline void sici_power_series(int sgn, std::complex<double> z,
						       std::complex<double> *s, std::complex<double> *c) {
	/* DLMF 6.6.5 and 6.6.6. If sgn = -1 computes si/ci, and if sgn = 1
	 * computes shi/chi.
	 */        
	std::complex<double> fac = z;
	*s = fac;
	*c = 0;
	std::complex<double> term1, term2;
	for (int n = 1; n < 100; n++) {
	    fac *= static_cast<double>(sgn)*z/(2.0*n);
	    term2 = fac/(2.0*n);
	    *c += term2;
	    fac *= z/(2.0*n + 1.0);
	    term1 = fac/(2.0*n + 1.0);
	    *s += term1;
	    constexpr double tol = std::numeric_limits<double>::epsilon();
	    if (std::abs(term1) < tol*std::abs(*s) && std::abs(term2) < tol*std::abs(*c)) {
		break;
	    }
	}
    }

}

    
XSF_HOST_DEVICE inline int sici(std::complex<double> z,
				    std::complex<double> *si, std::complex<double> *ci) {
    /* Compute sin/cos integrals at complex arguments. The algorithm
     * largely follows that of [1].
     */

    constexpr double EULER = xsf::cephes::detail::SCIPY_EULER;

    if (z == std::numeric_limits<double>::infinity()) {
        *si = M_PI_2;
        *ci = 0;
        return 0;
    }
    if (z == -std::numeric_limits<double>::infinity()) {
	*si = -M_PI_2;
        *ci = {0.0, M_PI};
        return 0;
    }

    if (std::abs(z) < 0.8) {
        // Use the series to avoid cancellation in si
	detail::sici_power_series(-1, z, si, ci);

        if (z == 0.0) {
            set_error("sici", SF_ERROR_DOMAIN, NULL);
            *ci = {-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()};
        } else {
            *ci += EULER + std::log(z);
	}
        return 0;
    }
    
    // DLMF 6.5.5/6.5.6 plus DLMF 6.4.4/6.4.6/6.4.7
    std::complex<double> jz = std::complex<double>(0.0, 1.0) * z;
    std::complex<double> term1 = expi(jz);
    std::complex<double> term2 = expi(-jz);
    *si = std::complex<double>(0.0, -0.5)*(term1 - term2);
    *ci = 0.5*(term1 + term2);
    if (z.real() == 0) {
        if (z.imag() > 0) {
            *ci += std::complex<double>(0.0, M_PI_2);
	} else if (z.imag() < 0) {
            *ci -= std::complex<double>(0.0, M_PI_2);
	}
    } else if (z.real() > 0) {
        *si -= M_PI_2;
    } else {
        *si += M_PI_2;
        if (z.imag() >= 0) {
            *ci += std::complex<double>(0.0, M_PI);
        } else {
            *ci -= std::complex<double>(0.0, M_PI);
	}
    }
    return 0;
}

XSF_HOST_DEVICE inline int sici(std::complex<float> z,
				    std::complex<float> *si_f, std::complex<float> *ci_f) {
    std::complex<double> si;
    std::complex<double> ci;
    int res = sici(z, &si, &ci);
    *si_f = si;
    *ci_f = ci;
    return res;
}

XSF_HOST_DEVICE inline int shichi(std::complex<double> z,
				       std::complex<double> *shi, std::complex<double> *chi) {
    /* Compute sinh/cosh integrals at complex arguments. The algorithm
     * largely follows that of [1].
     */
    constexpr double EULER = xsf::cephes::detail::SCIPY_EULER;
    if (z == std::numeric_limits<double>::infinity()) {
        *shi = std::numeric_limits<double>::infinity();
        *chi = std::numeric_limits<double>::infinity();
        return 0;
    }
    if (z == -std::numeric_limits<double>::infinity()) {
        *shi = -std::numeric_limits<double>::infinity();
        *chi = std::numeric_limits<double>::infinity();
        return 0;
    }
    if (std::abs(z) < 0.8) {
        // Use the series to avoid cancellation in shi
	detail::sici_power_series(1, z, shi, chi);
        if (z == 0.0) {
            set_error("shichi", SF_ERROR_DOMAIN, NULL);
            *chi = {-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()};
        } else {
            *chi += EULER + std::log(z);
	}
	return 0;
    }

    std::complex<double> term1 = expi(z);
    std::complex<double> term2 = expi(-z);
    *shi = 0.5*(term1 - term2);
    *chi = 0.5*(term1 + term2);
    if (z.imag() > 0) {
        *shi -= std::complex<double>(0.0, 0.5*M_PI);
        *chi += std::complex<double>(0.0, 0.5*M_PI);
    } else if (z.imag() < 0) {
        *shi += std::complex<double>(0.0, 0.5*M_PI);
        *chi -= std::complex<double>(0.0, 0.5*M_PI);
    } else if (z.real() < 0) {
        *chi += std::complex<double>(0.0, M_PI);
    }
    return 0;
}

XSF_HOST_DEVICE inline int shichi(std::complex<float> z,
				    std::complex<float> *shi_f, std::complex<float> *chi_f) {
    std::complex<double> shi;
    std::complex<double> chi;
    int res = shichi(z, &shi, &chi);
    *shi_f = shi;
    *chi_f = chi;
    return res;
}

XSF_HOST_DEVICE inline int sici(double x, double *si, double *ci) {
    return cephes::sici(x, si, ci);
}

XSF_HOST_DEVICE inline int shichi(double x, double *shi, double *chi) {
    return cephes::shichi(x, shi, chi);
}

XSF_HOST_DEVICE inline int sici(float x, float *si_f, float *ci_f) {
    double si;
    double ci;
    int res = cephes::sici(x, &si, &ci);
    *si_f = si;
    *ci_f = ci;
    return res;
}

XSF_HOST_DEVICE inline int shichi(float x, float *shi_f, float *chi_f) {
    double shi;
    double chi;
    int res = cephes::shichi(x, &shi, &chi);
    *shi_f = shi;
    *chi_f = chi;
    return res;
}
}
