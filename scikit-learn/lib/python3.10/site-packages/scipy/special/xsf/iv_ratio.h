// Numerically stable computation of iv(v+1, x) / iv(v, x)

#pragma once

#include "config.h"
#include "tools.h"
#include "error.h"
#include "cephes/dd_real.h"

namespace xsf {

/* Generates the "tail" of Perron's continued fraction for iv(v,x)/iv(v-1,x).
 *
 * The Perron continued fraction is studied in [1].  It is given by
 *
 *         iv(v, x)      x    -(2v+1)x   -(2v+3)x   -(2v+5)x
 *   R := --------- = ------ ---------- ---------- ---------- ...
 *        iv(v-1,x)   x+2v + 2(v+x)+1 + 2(v+x)+2 + 2(v+x)+3 +
 *
 * Given a suitable constant c, the continued fraction may be rearranged
 * into the following form to avoid premature floating point overflow:
 *
 *        xc                -(2vc+c)(xc) -(2vc+3c)(xc) -(2vc+5c)(xc)
 *   R = -----,  fc = 2vc + ------------ ------------- ------------- ...
 *       xc+fc              2(vc+xc)+c + 2(vc+xc)+2c + 2(vc+xc)+3c +
 *
 * This class generates the fractions of fc after 2vc.
 *
 * [1] Gautschi, W. and Slavik, J. (1978). "On the computation of modified
 *     Bessel function ratios." Mathematics of Computation, 32(143):865-875.
 */
template <class T>
struct IvRatioCFTailGenerator {

    XSF_HOST_DEVICE IvRatioCFTailGenerator(T vc, T xc, T c) noexcept {
        a0_ = -(2*vc-c)*xc;
        as_ = -2*c*xc;
        b0_ = 2*(vc+xc);
        bs_ = c;
        k_ = 0;
    }

    XSF_HOST_DEVICE std::pair<T, T> operator()() noexcept {
        using std::fma;
        ++k_;
        return {fma(static_cast<T>(k_), as_, a0_),
                fma(static_cast<T>(k_), bs_, b0_)};
    }

private:
    T a0_, as_;  // a[k] == a0 + as*k, k >= 1
    T b0_, bs_;  // b[k] == b0 + bs*k, k >= 1
    std::uint64_t k_; // current index
};

// Computes f(v, x) using Perron's continued fraction.
//
// T specifies the working type.  This allows the function to perform
// calculations in a higher precision, such as double-double, even if
// the return type is hardcoded to be double.
template <class T>
XSF_HOST_DEVICE inline std::pair<double, std::uint64_t>
_iv_ratio_cf(double v, double x, bool complement) {

    int e;
    std::frexp(std::fmax(v, x), &e);
    T c = T(std::ldexp(1, 2-e)); // rescaling multiplier
    T vc = v * c;
    T xc = x * c;

    IvRatioCFTailGenerator<T> cf(vc, xc, c);
    auto [fc, terms] = detail::series_eval_kahan(
        detail::continued_fraction_series(cf),
        T(std::numeric_limits<double>::epsilon()),
        1000,
        2*vc);

    T ret = (complement ? fc : xc) / (xc + fc);
    return {static_cast<double>(ret), terms};
}

XSF_HOST_DEVICE inline double iv_ratio(double v, double x) {

    if (std::isnan(v) || std::isnan(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (v < 0.5 || x < 0) {
        set_error("iv_ratio", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::isinf(v) && std::isinf(x)) {
        // There is not a unique limit as both v and x tends to infinity.
        set_error("iv_ratio", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (x == 0.0) {
        return x; // keep sign of x because iv_ratio is an odd function
    }
    if (std::isinf(v)) {
        return 0.0;
    }
    if (std::isinf(x)) {
        return 1.0;
    }

    auto [ret, terms] = _iv_ratio_cf<double>(v, x, false);
    if (terms == 0) { // failed to converge; should not happen
        set_error("iv_ratio", SF_ERROR_NO_RESULT, NULL);
        return std::numeric_limits<double>::quiet_NaN();
    }
    return ret;
}

XSF_HOST_DEVICE inline float iv_ratio(float v, float x) {
    return iv_ratio(static_cast<double>(v), static_cast<double>(x));
}

XSF_HOST_DEVICE inline double iv_ratio_c(double v, double x) {

    if (std::isnan(v) || std::isnan(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (v < 0.5 || x < 0) {
        set_error("iv_ratio_c", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::isinf(v) && std::isinf(x)) {
        // There is not a unique limit as both v and x tends to infinity.
        set_error("iv_ratio_c", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (x == 0.0) {
        return 1.0;
    }
    if (std::isinf(v)) {
        return 1.0;
    }
    if (std::isinf(x)) {
        return 0.0;
    }

    if (v >= 1) {
        // Numerical experiments show that evaluating the Perron c.f.
        // in double precision is sufficiently accurate if v >= 1.
        auto [ret, terms] = _iv_ratio_cf<double>(v, x, true);
        if (terms == 0) { // failed to converge; should not happen
            set_error("iv_ratio_c", SF_ERROR_NO_RESULT, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }
        return ret;
    } else if (v > 0.5) {
        // double-double arithmetic is needed for 0.5 < v < 1 to
        // achieve relative error on the scale of machine precision.
        using cephes::detail::double_double;
        auto [ret, terms] = _iv_ratio_cf<double_double>(v, x, true);
        if (terms == 0) { // failed to converge; should not happen
            set_error("iv_ratio_c", SF_ERROR_NO_RESULT, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }
        return ret;
    } else {
        // The previous branch (v > 0.5) also works for v == 0.5, but
        // the closed-form formula "1 - tanh(x)" is more efficient.
        double t = std::exp(-2*x);
        return (2 * t) / (1 + t);
    }
}

XSF_HOST_DEVICE inline float iv_ratio_c(float v, float x) {
    return iv_ratio_c(static_cast<double>(v), static_cast<double>(x));
}

} // namespace xsf
