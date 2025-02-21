/* Implementation of Gauss's hypergeometric function for complex values.
 *
 * This implementation is based on the Fortran implementation by Shanjie Zhang and
 * Jianming Jin included in specfun.f [1]_.  Computation of Gauss's hypergeometric
 * function involves handling a patchwork of special cases. By default the Zhang and
 * Jin implementation has been followed as closely as possible except for situations where
 * an improvement was obvious. We've attempted to document the reasons behind decisions
 * made by Zhang and Jin and to document the reasons for deviating from their implementation
 * when this has been done. References to the NIST Digital Library of Mathematical
 * Functions [2]_ have been added where they are appropriate. The review paper by
 * Pearson et al [3]_ is an excellent resource for best practices for numerical
 * computation of hypergeometric functions. We have followed this review paper
 * when making improvements to and correcting defects in Zhang and Jin's
 * implementation. When Pearson et al propose several competing alternatives for a
 * given case, we've used our best judgment to decide on the method to use.
 *
 * Author: Albert Steppi
 *
 * Distributed under the same license as Scipy.
 *
 * References
 * ----------
 * .. [1] S. Zhang and J.M. Jin, "Computation of Special Functions", Wiley 1996
 * .. [2] NIST Digital Library of Mathematical Functions. http://dlmf.nist.gov/,
 *        Release 1.1.1 of 2021-03-15. F. W. J. Olver, A. B. Olde Daalhuis,
 *        D. W. Lozier, B. I. Schneider, R. F. Boisvert, C. W. Clark, B. R. Miller,
 *        B. V. Saunders, H. S. Cohl, and M. A. McClain, eds.
 * .. [3] Pearson, J.W., Olver, S. & Porter, M.A.
 *        "Numerical methods for the computation of the confluent and Gauss
 *        hypergeometric functions."
 *        Numer Algor 74, 821-866 (2017). https://doi.org/10.1007/s11075-016-0173-0
 * .. [4] Raimundas Vidunas, "Degenerate Gauss Hypergeometric Functions",
 *        Kyushu Journal of Mathematics, 2007, Volume 61, Issue 1, Pages 109-135,
 * .. [5] López, J.L., Temme, N.M. New series expansions of the Gauss hypergeometric
 *        function. Adv Comput Math 39, 349-365 (2013).
 *        https://doi.org/10.1007/s10444-012-9283-y
 * """
 */

#pragma once

#include "config.h"
#include "error.h"
#include "tools.h"

#include "binom.h"
#include "cephes/gamma.h"
#include "cephes/lanczos.h"
#include "cephes/poch.h"
#include "cephes/hyp2f1.h"
#include "digamma.h"

namespace xsf {
namespace detail {
    constexpr double hyp2f1_EPS = 1e-15;
    /* The original implementation in SciPy from Zhang and Jin used 1500 for the
     * maximum number of series iterations in some cases and 500 in others.
     * Through the empirical results on the test cases in
     * scipy/special/_precompute/hyp2f1_data.py, it was determined that these values
     * can lead to early termination of series which would have eventually converged
     * at a reasonable level of accuracy. We've bumped the iteration limit to 3000,
     * and may adjust it again based on further analysis. */
    constexpr std::uint64_t hyp2f1_MAXITER = 3000;

    XSF_HOST_DEVICE inline double four_gammas_lanczos(double u, double v, double w, double x) {
        /* Compute ratio of gamma functions using lanczos approximation.
         *
         * Computes gamma(u)*gamma(v)/(gamma(w)*gamma(x))
         *
         * It is assumed that x = u + v - w, but it is left to the user to
         * ensure this.
         *
         * The lanczos approximation takes the form
         *
         * gamma(x) = factor(x) * lanczos_sum_expg_scaled(x)
         *
         * where factor(x) = ((x + lanczos_g - 0.5)/e)**(x - 0.5).
         *
         * The formula above is only valid for x >= 0.5, but can be extended to
         * x < 0.5 with the reflection principle.
         *
         * Using the lanczos approximation when computing this ratio of gamma functions
         * allows factors to be combined analytically to avoid underflow and overflow
         * and produce a more accurate result. The condition x = u + v - w makes it
         * possible to cancel the factors in the expression
         *
         * factor(u) * factor(v) / (factor(w) * factor(x))
         *
         * by taking one factor and absorbing it into the others. Currently, this
         * implementation takes the factor corresponding to the argument with largest
         * absolute value and absorbs it into the others.
         *
         * Since this is only called internally by four_gammas. It is assumed that
         * |u| >= |v| and |w| >= |x|.
         */

        /* The below implementation may incorrectly return finite results
         * at poles of the gamma function. Handle these cases explicitly. */
        if ((u == std::trunc(u) && u <= 0) || (v == std::trunc(v) && v <= 0)) {
            /* Return nan if numerator has pole. Diverges to +- infinity
             * depending on direction so value is undefined. */
            return std::numeric_limits<double>::quiet_NaN();
        }
        if ((w == std::trunc(w) && w <= 0) || (x == std::trunc(x) && x <= 0)) {
            // Return 0 if denominator has pole but not numerator.
            return 0.0;
        }

        double result = 1.0;
        double ugh, vgh, wgh, xgh, u_prime, v_prime, w_prime, x_prime;

        if (u >= 0.5) {
            result *= cephes::lanczos_sum_expg_scaled(u);
            ugh = u + cephes::lanczos_g - 0.5;
            u_prime = u;
        } else {
            result /= cephes::lanczos_sum_expg_scaled(1 - u) * std::sin(M_PI * u) * M_1_PI;
            ugh = 0.5 - u + cephes::lanczos_g;
            u_prime = 1 - u;
        }

        if (v >= 0.5) {
            result *= cephes::lanczos_sum_expg_scaled(v);
            vgh = v + cephes::lanczos_g - 0.5;
            v_prime = v;
        } else {
            result /= cephes::lanczos_sum_expg_scaled(1 - v) * std::sin(M_PI * v) * M_1_PI;
            vgh = 0.5 - v + cephes::lanczos_g;
            v_prime = 1 - v;
        }

        if (w >= 0.5) {
            result /= cephes::lanczos_sum_expg_scaled(w);
            wgh = w + cephes::lanczos_g - 0.5;
            w_prime = w;
        } else {
            result *= cephes::lanczos_sum_expg_scaled(1 - w) * std::sin(M_PI * w) * M_1_PI;
            wgh = 0.5 - w + cephes::lanczos_g;
            w_prime = 1 - w;
        }

        if (x >= 0.5) {
            result /= cephes::lanczos_sum_expg_scaled(x);
            xgh = x + cephes::lanczos_g - 0.5;
            x_prime = x;
        } else {
            result *= cephes::lanczos_sum_expg_scaled(1 - x) * std::sin(M_PI * x) * M_1_PI;
            xgh = 0.5 - x + cephes::lanczos_g;
            x_prime = 1 - x;
        }

        if (std::abs(u) >= std::abs(w)) {
            // u has greatest absolute value. Absorb ugh into the others.
            if (std::abs((v_prime - u_prime) * (v - 0.5)) < 100 * ugh and v > 100) {
                /* Special case where base is close to 1. Condition taken from
                 * Boost's beta function implementation. */
                result *= std::exp((v - 0.5) * std::log1p((v_prime - u_prime) / ugh));
            } else {
                result *= std::pow(vgh / ugh, v - 0.5);
            }

            if (std::abs((u_prime - w_prime) * (w - 0.5)) < 100 * wgh and u > 100) {
                result *= std::exp((w - 0.5) * std::log1p((u_prime - w_prime) / wgh));
            } else {
                result *= std::pow(ugh / wgh, w - 0.5);
            }

            if (std::abs((u_prime - x_prime) * (x - 0.5)) < 100 * xgh and u > 100) {
                result *= std::exp((x - 0.5) * std::log1p((u_prime - x_prime) / xgh));
            } else {
                result *= std::pow(ugh / xgh, x - 0.5);
            }
        } else {
            // w has greatest absolute value. Absorb wgh into the others.
            if (std::abs((u_prime - w_prime) * (u - 0.5)) < 100 * wgh and u > 100) {
                result *= std::exp((u - 0.5) * std::log1p((u_prime - w_prime) / wgh));
            } else {
                result *= pow(ugh / wgh, u - 0.5);
            }
            if (std::abs((v_prime - w_prime) * (v - 0.5)) < 100 * wgh and v > 100) {
                result *= std::exp((v - 0.5) * std::log1p((v_prime - w_prime) / wgh));
            } else {
                result *= std::pow(vgh / wgh, v - 0.5);
            }
            if (std::abs((w_prime - x_prime) * (x - 0.5)) < 100 * xgh and x > 100) {
                result *= std::exp((x - 0.5) * std::log1p((w_prime - x_prime) / xgh));
            } else {
                result *= std::pow(wgh / xgh, x - 0.5);
            }
        }
        // This exhausts all cases because we assume |u| >= |v| and |w| >= |x|.

        return result;
    }

    XSF_HOST_DEVICE inline double four_gammas(double u, double v, double w, double x) {
        double result;

        // Without loss of generality, ensure |u| >= |v| and |w| >= |x|.
        if (std::abs(v) > std::abs(u)) {
            std::swap(u, v);
        }
        if (std::abs(x) > std::abs(w)) {
            std::swap(x, w);
        }
        /* Direct ratio tends to be more accurate for arguments in this range. Range
         * chosen empirically based on the relevant benchmarks in
         * scipy/special/_precompute/hyp2f1_data.py */
        if (std::abs(u) <= 100 && std::abs(v) <= 100 && std::abs(w) <= 100 && std::abs(x) <= 100) {
            result = cephes::Gamma(u) * cephes::Gamma(v) * (cephes::rgamma(w) * cephes::rgamma(x));
            if (std::isfinite(result) && result != 0.0) {
                return result;
            }
        }
        result = four_gammas_lanczos(u, v, w, x);
        if (std::isfinite(result) && result != 0.0) {
            return result;
        }
        // If overflow or underflow, try again with logs.
        result = std::exp(cephes::lgam(v) - cephes::lgam(x) + cephes::lgam(u) - cephes::lgam(w));
        result *= cephes::gammasgn(u) * cephes::gammasgn(w) * cephes::gammasgn(v) * cephes::gammasgn(x);
        return result;
    }

    class HypergeometricSeriesGenerator {
        /* Maclaurin series for hyp2f1.
         *
         * Series is convergent for |z| < 1 but is only practical for numerical
         * computation when |z| < 0.9.
         */
      public:
        XSF_HOST_DEVICE HypergeometricSeriesGenerator(double a, double b, double c, std::complex<double> z)
            : a_(a), b_(b), c_(c), z_(z), term_(1.0), k_(0) {}

        XSF_HOST_DEVICE std::complex<double> operator()() {
            std::complex<double> output = term_;
            term_ = term_ * (a_ + k_) * (b_ + k_) / ((k_ + 1) * (c_ + k_)) * z_;
            ++k_;
            return output;
        }

      private:
        double a_, b_, c_;
        std::complex<double> z_, term_;
        std::uint64_t k_;
    };

    class Hyp2f1Transform1Generator {
        /* 1 -z transformation of standard series.*/
      public:
        XSF_HOST_DEVICE Hyp2f1Transform1Generator(double a, double b, double c, std::complex<double> z)
            : factor1_(four_gammas(c, c - a - b, c - a, c - b)),
              factor2_(four_gammas(c, a + b - c, a, b) * std::pow(1.0 - z, c - a - b)),
              generator1_(HypergeometricSeriesGenerator(a, b, a + b - c + 1, 1.0 - z)),
              generator2_(HypergeometricSeriesGenerator(c - a, c - b, c - a - b + 1, 1.0 - z)) {}

        XSF_HOST_DEVICE std::complex<double> operator()() {
            return factor1_ * generator1_() + factor2_ * generator2_();
        }

      private:
        std::complex<double> factor1_, factor2_;
        HypergeometricSeriesGenerator generator1_, generator2_;
    };

    class Hyp2f1Transform1LimitSeriesGenerator {
        /* 1 - z transform in limit as c - a - b approaches an integer m. */
      public:
        XSF_HOST_DEVICE Hyp2f1Transform1LimitSeriesGenerator(double a, double b, double m, std::complex<double> z)
            : d1_(xsf::digamma(a)), d2_(xsf::digamma(b)), d3_(xsf::digamma(1 + m)),
              d4_(xsf::digamma(1.0)), a_(a), b_(b), m_(m), z_(z), log_1_z_(std::log(1.0 - z)),
              factor_(cephes::rgamma(m + 1)), k_(0) {}

        XSF_HOST_DEVICE std::complex<double> operator()() {
            std::complex<double> term_ = (d1_ + d2_ - d3_ - d4_ + log_1_z_) * factor_;
            // Use digamma(x + 1) = digamma(x) + 1/x
            d1_ += 1 / (a_ + k_);       // d1 = digamma(a + k)
            d2_ += 1 / (b_ + k_);       // d2 = digamma(b + k)
            d3_ += 1 / (1.0 + m_ + k_); // d3 = digamma(1 + m + k)
            d4_ += 1 / (1.0 + k_);      // d4 = digamma(1 + k)
            factor_ *= (a_ + k_) * (b_ + k_) / ((k_ + 1.0) * (m_ + k_ + 1)) * (1.0 - z_);
            ++k_;
            return term_;
        }

      private:
        double d1_, d2_, d3_, d4_, a_, b_, m_;
        std::complex<double> z_, log_1_z_, factor_;
        int k_;
    };

    class Hyp2f1Transform2Generator {
        /* 1/z transformation of standard series.*/
      public:
        XSF_HOST_DEVICE Hyp2f1Transform2Generator(double a, double b, double c, std::complex<double> z)
            : factor1_(four_gammas(c, b - a, b, c - a) * std::pow(-z, -a)),
              factor2_(four_gammas(c, a - b, a, c - b) * std::pow(-z, -b)),
              generator1_(HypergeometricSeriesGenerator(a, a - c + 1, a - b + 1, 1.0 / z)),
              generator2_(HypergeometricSeriesGenerator(b, b - c + 1, b - a + 1, 1.0 / z)) {}

        XSF_HOST_DEVICE std::complex<double> operator()() {
            return factor1_ * generator1_() + factor2_ * generator2_();
        }

      private:
        std::complex<double> factor1_, factor2_;
        HypergeometricSeriesGenerator generator1_, generator2_;
    };

    class Hyp2f1Transform2LimitSeriesGenerator {
        /* 1/z transform in limit as a - b approaches a non-negative integer m. (Can swap a and b to
         * handle the m a negative integer case. */
      public:
        XSF_HOST_DEVICE Hyp2f1Transform2LimitSeriesGenerator(double a, double b, double c, double m,
                                                                 std::complex<double> z)
            : d1_(xsf::digamma(1.0)), d2_(xsf::digamma(1 + m)), d3_(xsf::digamma(a)),
              d4_(xsf::digamma(c - a)), a_(a), b_(b), c_(c), m_(m), z_(z), log_neg_z_(std::log(-z)),
              factor_(xsf::cephes::poch(b, m) * xsf::cephes::poch(1 - c + b, m) *
                      xsf::cephes::rgamma(m + 1)),
              k_(0) {}

        XSF_HOST_DEVICE std::complex<double> operator()() {
            std::complex<double> term = (d1_ + d2_ - d3_ - d4_ + log_neg_z_) * factor_;
            // Use digamma(x + 1) = digamma(x) + 1/x
            d1_ += 1 / (1.0 + k_);         // d1 = digamma(1 + k)
            d2_ += 1 / (1.0 + m_ + k_);    // d2 = digamma(1 + m + k)
            d3_ += 1 / (a_ + k_);          // d3 = digamma(a + k)
            d4_ -= 1 / (c_ - a_ - k_ - 1); // d4 = digamma(c - a - k)
            factor_ *= (b_ + m_ + k_) * (1 - c_ + b_ + m_ + k_) / ((k_ + 1) * (m_ + k_ + 1)) / z_;
            ++k_;
            return term;
        }

      private:
        double d1_, d2_, d3_, d4_, a_, b_, c_, m_;
        std::complex<double> z_, log_neg_z_, factor_;
        std::uint64_t k_;
    };

    class Hyp2f1Transform2LimitSeriesCminusAIntGenerator {
        /* 1/z transform in limit as a - b approaches a non-negative integer m, and c - a approaches
         * a positive integer n. */
      public:
        XSF_HOST_DEVICE Hyp2f1Transform2LimitSeriesCminusAIntGenerator(double a, double b, double c, double m,
                                                                           double n, std::complex<double> z)
            : d1_(xsf::digamma(1.0)), d2_(xsf::digamma(1 + m)), d3_(xsf::digamma(a)),
              d4_(xsf::digamma(n)), a_(a), b_(b), c_(c), m_(m), n_(n), z_(z), log_neg_z_(std::log(-z)),
              factor_(xsf::cephes::poch(b, m) * xsf::cephes::poch(1 - c + b, m) *
                      xsf::cephes::rgamma(m + 1)),
              k_(0) {}

        XSF_HOST_DEVICE std::complex<double> operator()() {
            std::complex<double> term;
            if (k_ < n_) {
                term = (d1_ + d2_ - d3_ - d4_ + log_neg_z_) * factor_;
                // Use digamma(x + 1) = digamma(x) + 1/x
                d1_ += 1 / (1.0 + k_);    // d1 = digamma(1 + k)
                d2_ += 1 / (1 + m_ + k_); // d2 = digamma(1 + m + k)
                d3_ += 1 / (a_ + k_);     // d3 = digamma(a + k)
                d4_ -= 1 / (n_ - k_ - 1); // d4 = digamma(c - a - k)
                factor_ *= (b_ + m_ + k_) * (1 - c_ + b_ + m_ + k_) / ((k_ + 1) * (m_ + k_ + 1)) / z_;
                ++k_;
                return term;
            }
            if (k_ == n_) {
                /* When c - a approaches a positive integer and k_ >= c - a = n then
                 * poch(1 - c + b + m + k) = poch(1 - c + a + k) = approaches zero and
                 * digamma(c - a - k) approaches a pole. However we can use the limit
                 * digamma(-n + epsilon) / gamma(-n + epsilon) -> (-1)**(n + 1) * (n+1)! as epsilon -> 0
                 * to continue the series.
                 *
                 * poch(1 - c + b, m + k) = gamma(1 - c + b + m + k)/gamma(1 - c + b)
                 *
                 * If a - b is an integer and c - a is an integer, then a and b must both be integers, so assume
                 * a and b are integers and take the limit as c approaches an integer.
                 *
                 * gamma(1 - c + epsilon + a + k)/gamma(1 - c - epsilon + b) =
                 * (gamma(c + epsilon - b) / gamma(c + epsilon - a - k)) *
                 * (sin(pi * (c + epsilon - b)) / sin(pi * (c + epsilon - a - k))) (reflection principle)
                 *
                 * In the limit as epsilon goes to zero, the ratio of sines will approach
                 * (-1)**(a - b + k) = (-1)**(m + k)
                 *
                 * We may then replace
                 *
                 * poch(1 - c - epsilon + b, m + k)*digamma(c + epsilon - a - k)
                 *
                 * with
                 *
                 * (-1)**(a - b + k)*gamma(c + epsilon - b) * digamma(c + epsilon - a - k) / gamma(c + epsilon - a - k)
                 *
                 * and taking the limit epsilon -> 0 gives
                 *
                 * (-1)**(a - b + k) * gamma(c - b) * (-1)**(k + a - c + 1)(k + a - c)!
                 * = (-1)**(c - b - 1)*Gamma(k + a - c + 1)
                 */
                factor_ = std::pow(-1, m_ + n_) * xsf::binom(c_ - 1, b_ - 1) *
                          xsf::cephes::poch(c_ - a_ + 1, m_ - 1) / std::pow(z_, static_cast<double>(k_));
            }
            term = factor_;
            factor_ *= (b_ + m_ + k_) * (k_ + a_ - c_ + 1) / ((k_ + 1) * (m_ + k_ + 1)) / z_;
            ++k_;
            return term;
        }

      private:
        double d1_, d2_, d3_, d4_, a_, b_, c_, m_, n_;
        std::complex<double> z_, log_neg_z_, factor_;
        std::uint64_t k_;
    };

    class Hyp2f1Transform2LimitFinitePartGenerator {
        /* Initial finite sum in limit as a - b approaches a non-negative integer m. The limiting series
         * for the 1 - z transform also has an initial finite sum, but it is a standard hypergeometric
         * series. */
      public:
        XSF_HOST_DEVICE Hyp2f1Transform2LimitFinitePartGenerator(double b, double c, double m,
                                                                     std::complex<double> z)
            : b_(b), c_(c), m_(m), z_(z), term_(cephes::Gamma(m) * cephes::rgamma(c - b)), k_(0) {}

        XSF_HOST_DEVICE std::complex<double> operator()() {
            std::complex<double> output = term_;
            term_ = term_ * (b_ + k_) * (c_ - b_ - k_ - 1) / ((k_ + 1) * (m_ - k_ - 1)) / z_;
            ++k_;
            return output;
        }

      private:
        double b_, c_, m_;
        std::complex<double> z_, term_;
        std::uint64_t k_;
    };

    class LopezTemmeSeriesGenerator {
        /* Lopez-Temme Series for Gaussian hypergeometric function [4].
         *
         * Converges for all z with real(z) < 1, including in the regions surrounding
         * the points exp(+- i*pi/3) that are not covered by any of the standard
         * transformations.
         */
      public:
        XSF_HOST_DEVICE LopezTemmeSeriesGenerator(double a, double b, double c, std::complex<double> z)
            : n_(0), a_(a), b_(b), c_(c), phi_previous_(1.0), phi_(1 - 2 * b / c), z_(z), Z_(a * z / (z - 2.0)) {}

        XSF_HOST_DEVICE std::complex<double> operator()() {
            if (n_ == 0) {
                ++n_;
                return 1.0;
            }
            if (n_ > 1) { // Update phi and Z for n>=2
                double new_phi = ((n_ - 1) * phi_previous_ - (2.0 * b_ - c_) * phi_) / (c_ + (n_ - 1));
                phi_previous_ = phi_;
                phi_ = new_phi;
                Z_ = Z_ * z_ / (z_ - 2.0) * ((a_ + (n_ - 1)) / n_);
            }
            ++n_;
            return Z_ * phi_;
        }

      private:
        std::uint64_t n_;
        double a_, b_, c_, phi_previous_, phi_;
        std::complex<double> z_, Z_;
    };

    XSF_HOST_DEVICE std::complex<double> hyp2f1_transform1_limiting_case(double a, double b, double c, double m,
                                                                             std::complex<double> z) {
        /* 1 - z transform in limiting case where c - a - b approaches an integer m. */
        std::complex<double> result = 0.0;
        if (m >= 0) {
            if (m != 0) {
                auto series_generator = HypergeometricSeriesGenerator(a, b, 1 - m, 1.0 - z);
                result += four_gammas(m, c, a + m, b + m) * series_eval_fixed_length(series_generator,
                                                                                     std::complex<double>{0.0, 0.0},
                                                                                     static_cast<std::uint64_t>(m));
            }
            std::complex<double> prefactor = std::pow(-1.0, m + 1) * xsf::cephes::Gamma(c) /
                                             (xsf::cephes::Gamma(a) * xsf::cephes::Gamma(b)) *
                                             std::pow(1.0 - z, m);
            auto series_generator = Hyp2f1Transform1LimitSeriesGenerator(a + m, b + m, m, z);
            result += prefactor * series_eval(series_generator, std::complex<double>{0.0, 0.0}, hyp2f1_EPS,
                                              hyp2f1_MAXITER, "hyp2f1");
            return result;
        } else {
            result = four_gammas(-m, c, a, b) * std::pow(1.0 - z, m);
            auto series_generator1 = HypergeometricSeriesGenerator(a + m, b + m, 1 + m, 1.0 - z);
            result *= series_eval_fixed_length(series_generator1, std::complex<double>{0.0, 0.0},
                                               static_cast<std::uint64_t>(-m));
            double prefactor = std::pow(-1.0, m + 1) * xsf::cephes::Gamma(c) *
                               (xsf::cephes::rgamma(a + m) * xsf::cephes::rgamma(b + m));
            auto series_generator2 = Hyp2f1Transform1LimitSeriesGenerator(a, b, -m, z);
            result += prefactor * series_eval(series_generator2, std::complex<double>{0.0, 0.0}, hyp2f1_EPS,
                                              hyp2f1_MAXITER, "hyp2f1");
            return result;
        }
    }

    XSF_HOST_DEVICE std::complex<double> hyp2f1_transform2_limiting_case(double a, double b, double c, double m,
                                                                             std::complex<double> z) {
        /* 1 / z transform in limiting case where a - b approaches a non-negative integer m. Negative integer case
         * can be handled by swapping a and b. */
        auto series_generator1 = Hyp2f1Transform2LimitFinitePartGenerator(b, c, m, z);
        std::complex<double> result = cephes::Gamma(c) * cephes::rgamma(a) * std::pow(-z, -b);
        result *=
            series_eval_fixed_length(series_generator1, std::complex<double>{0.0, 0.0}, static_cast<std::uint64_t>(m));
        std::complex<double> prefactor = cephes::Gamma(c) * (cephes::rgamma(a) * cephes::rgamma(c - b) * std::pow(-z, -a));
        double n = c - a;
        if (abs(n - std::round(n)) < hyp2f1_EPS) {
            auto series_generator2 = Hyp2f1Transform2LimitSeriesCminusAIntGenerator(a, b, c, m, n, z);
            result += prefactor * series_eval(series_generator2, std::complex<double>{0.0, 0.0}, hyp2f1_EPS,
                                              hyp2f1_MAXITER, "hyp2f1");
            return result;
        }
        auto series_generator2 = Hyp2f1Transform2LimitSeriesGenerator(a, b, c, m, z);
        result += prefactor *
                  series_eval(series_generator2, std::complex<double>{0.0, 0.0}, hyp2f1_EPS, hyp2f1_MAXITER, "hyp2f1");
        return result;
    }

} // namespace detail

XSF_HOST_DEVICE inline std::complex<double> hyp2f1(double a, double b, double c, std::complex<double> z) {
    /* Special Cases
     * -----------------------------------------------------------------------
     * Takes constant value 1 when a = 0 or b = 0, even if c is a non-positive
     * integer. This follows mpmath. */
    if (a == 0 || b == 0) {
        return 1.0;
    }
    double z_abs = std::abs(z);
    // Equals 1 when z i 0, unless c is 0.
    if (z_abs == 0) {
        if (c != 0) {
            return 1.0;
        } else {
            // Returning real part NAN and imaginary part 0 follows mpmath.
            return std::complex<double>{std::numeric_limits<double>::quiet_NaN(), 0};
        }
    }
    bool a_neg_int = a == std::trunc(a) && a < 0;
    bool b_neg_int = b == std::trunc(b) && b < 0;
    bool c_non_pos_int = c == std::trunc(c) and c <= 0;
    /* Diverges when c is a non-positive integer unless a is an integer with
     * c <= a <= 0 or b is an integer with c <= b <= 0, (or z equals 0 with
     * c != 0) Cases z = 0, a = 0, or b = 0 have already been handled. We follow
     * mpmath in handling the degenerate cases where any of a, b, c are
     * non-positive integers. See [3] for a treatment of degenerate cases. */
    if (c_non_pos_int && !((a_neg_int && c <= a && a < 0) || (b_neg_int && c <= b && b < 0))) {
        return std::complex<double>{std::numeric_limits<double>::infinity(), 0};
    }
    /* Reduces to a polynomial when a or b is a negative integer.
     * If a and b are both negative integers, we take care to terminate
     * the series at a or b of smaller magnitude. This is to ensure proper
     * handling of situations like a < c < b <= 0, a, b, c all non-positive
     * integers, where terminating at a would lead to a term of the form 0 / 0. */
    std::uint64_t max_degree;
    if (a_neg_int || b_neg_int) {
        if (a_neg_int && b_neg_int) {
            max_degree = a > b ? std::abs(a) : std::abs(b);
        } else if (a_neg_int) {
            max_degree = std::abs(a);
        } else {
            max_degree = std::abs(b);
        }
        if (max_degree <= UINT64_MAX) {
            auto series_generator = detail::HypergeometricSeriesGenerator(a, b, c, z);
            return detail::series_eval_fixed_length(series_generator, std::complex<double>{0.0, 0.0}, max_degree + 1);
        } else {
            set_error("hyp2f1", SF_ERROR_NO_RESULT, NULL);
            return std::complex<double>{std::numeric_limits<double>::quiet_NaN(),
                                        std::numeric_limits<double>::quiet_NaN()};
        }
    }
    // Kummer's Theorem for z = -1; c = 1 + a - b (DLMF 15.4.26)
    if (std::abs(z + 1.0) < detail::hyp2f1_EPS && std::abs(1 + a - b - c) < detail::hyp2f1_EPS && !c_non_pos_int) {
        return detail::four_gammas(a - b + 1, 0.5 * a + 1, a + 1, 0.5 * a - b + 1);
    }
    std::complex<double> result;
    bool c_minus_a_neg_int = c - a == std::trunc(c - a) && c - a < 0;
    bool c_minus_b_neg_int = c - b == std::trunc(c - b) && c - b < 0;
    /* If one of c - a or c - b is a negative integer, reduces to evaluating
     * a polynomial through an Euler hypergeometric transformation.
     * (DLMF 15.8.1) */
    if (c_minus_a_neg_int || c_minus_b_neg_int) {
        max_degree = c_minus_b_neg_int ? std::abs(c - b) : std::abs(c - a);
        if (max_degree <= UINT64_MAX) {
            result = std::pow(1.0 - z, c - a - b);
            auto series_generator = detail::HypergeometricSeriesGenerator(c - a, c - b, c, z);
            result *=
                detail::series_eval_fixed_length(series_generator, std::complex<double>{0.0, 0.0}, max_degree + 2);
            return result;
        } else {
            set_error("hyp2f1", SF_ERROR_NO_RESULT, NULL);
            return std::complex<double>{std::numeric_limits<double>::quiet_NaN(),
                                        std::numeric_limits<double>::quiet_NaN()};
        }
    }
    /* Diverges as real(z) -> 1 when c <= a + b.
     * Todo: Actually check for overflow instead of using a fixed tolerance for
     * all parameter combinations like in the Fortran original. */
    if (std::abs(1 - z.real()) < detail::hyp2f1_EPS && z.imag() == 0 && c - a - b <= 0 && !c_non_pos_int) {
        return std::complex<double>{std::numeric_limits<double>::infinity(), 0};
    }
    // Gauss's Summation Theorem for z = 1; c - a - b > 0 (DLMF 15.4.20).
    if (z == 1.0 && c - a - b > 0 && !c_non_pos_int) {
        return detail::four_gammas(c, c - a - b, c - a, c - b);
    }
    /* |z| < 0, z.real() >= 0. Use the Maclaurin Series.
     * -----------------------------------------------------------------------
     * Apply Euler Hypergeometric Transformation (DLMF 15.8.1) to reduce
     * size of a and b if possible. We follow Zhang and Jin's
     * implementation [1] although there is very likely a better heuristic
     * to determine when this transformation should be applied. As it
     * stands, this hurts precision in some cases. */
    if (z_abs < 0.9 && z.real() >= 0) {
        if (c - a < a && c - b < b) {
            result = std::pow(1.0 - z, c - a - b);
            auto series_generator = detail::HypergeometricSeriesGenerator(c - a, c - b, c, z);
            result *= detail::series_eval(series_generator, std::complex<double>{0.0, 0.0}, detail::hyp2f1_EPS,
                                          detail::hyp2f1_MAXITER, "hyp2f1");
            return result;
        }
        auto series_generator = detail::HypergeometricSeriesGenerator(a, b, c, z);
        return detail::series_eval(series_generator, std::complex<double>{0.0, 0.0}, detail::hyp2f1_EPS,
                                   detail::hyp2f1_MAXITER, "hyp2f1");
    }
    /* Points near exp(iπ/3), exp(-iπ/3) not handled by any of the standard
     * transformations. Use series of López and Temme [5]. These regions
     * were not correctly handled by Zhang and Jin's implementation.
     * -------------------------------------------------------------------------*/
    if (0.9 <= z_abs && z_abs < 1.1 && std::abs(1.0 - z) >= 0.9 && z.real() >= 0) {
        /* This condition for applying Euler Transformation (DLMF 15.8.1)
         * was determined empirically to work better for this case than that
         * used in Zhang and Jin's implementation for |z| < 0.9,
         *  real(z) >= 0. */
        if ((c - a <= a && c - b < b) || (c - a < a && c - b <= b)) {
            auto series_generator = detail::LopezTemmeSeriesGenerator(c - a, c - b, c, z);
            result = std::pow(1.0 - 0.5 * z, a - c); // Lopez-Temme prefactor
            result *= detail::series_eval(series_generator, std::complex<double>{0.0, 0.0}, detail::hyp2f1_EPS,
                                          detail::hyp2f1_MAXITER, "hyp2f1");
            return std::pow(1.0 - z, c - a - b) * result; // Euler transform prefactor.
        }
        auto series_generator = detail::LopezTemmeSeriesGenerator(a, b, c, z);
        result = detail::series_eval(series_generator, std::complex<double>{0.0, 0.0}, detail::hyp2f1_EPS,
                                     detail::hyp2f1_MAXITER, "hyp2f1");
        return std::pow(1.0 - 0.5 * z, -a) * result; // Lopez-Temme prefactor.
    }
    /* z/(z - 1) transformation (DLMF 15.8.1). Avoids cancellation issues that
     * occur with Maclaurin series for real(z) < 0.
     * -------------------------------------------------------------------------*/
    if (z_abs < 1.1 && z.real() < 0) {
        if (0 < b && b < a && a < c) {
            std::swap(a, b);
        }
        auto series_generator = detail::HypergeometricSeriesGenerator(a, c - b, c, z / (z - 1.0));
        return std::pow(1.0 - z, -a) * detail::series_eval(series_generator, std::complex<double>{0.0, 0.0},
                                                           detail::hyp2f1_EPS, detail::hyp2f1_MAXITER, "hyp2f1");
    }
    /* 1 - z transformation (DLMF 15.8.4). */
    if (0.9 <= z_abs && z_abs < 1.1) {
        if (std::abs(c - a - b - std::round(c - a - b)) < detail::hyp2f1_EPS) {
            // Removable singularity when c - a - b is an integer. Need to use limiting formula.
            double m = std::round(c - a - b);
            return detail::hyp2f1_transform1_limiting_case(a, b, c, m, z);
        }
        auto series_generator = detail::Hyp2f1Transform1Generator(a, b, c, z);
        return detail::series_eval(series_generator, std::complex<double>{0.0, 0.0}, detail::hyp2f1_EPS,
                                   detail::hyp2f1_MAXITER, "hyp2f1");
    }
    /* 1/z transformation (DLMF 15.8.2). */
    if (std::abs(a - b - std::round(a - b)) < detail::hyp2f1_EPS) {
        if (b > a) {
            std::swap(a, b);
        }
        double m = std::round(a - b);
        return detail::hyp2f1_transform2_limiting_case(a, b, c, m, z);
    }
    auto series_generator = detail::Hyp2f1Transform2Generator(a, b, c, z);
    return detail::series_eval(series_generator, std::complex<double>{0.0, 0.0}, detail::hyp2f1_EPS,
                               detail::hyp2f1_MAXITER, "hyp2f1");
}

XSF_HOST_DEVICE inline std::complex<float> hyp2f1(float a, float b, float c, std::complex<float> x) {
    return static_cast<std::complex<float>>(hyp2f1(static_cast<double>(a), static_cast<double>(b),
                                                   static_cast<double>(c), static_cast<std::complex<double>>(x)));
}

XSF_HOST_DEVICE inline double hyp2f1(double a, double b, double c, double x) { return cephes::hyp2f1(a, b, c, x); }

XSF_HOST_DEVICE inline float hyp2f1(float a, float b, float c, float x) {
    return hyp2f1(static_cast<double>(a), static_cast<double>(b), static_cast<double>(c), static_cast<double>(x));
}

} // namespace xsf
