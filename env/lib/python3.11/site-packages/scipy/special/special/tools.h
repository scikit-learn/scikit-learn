/* Building blocks for implementing special functions */

#pragma once

#include "config.h"
#include "error.h"

namespace special {
namespace detail {

    /* Result type of a "generator", a callable object that produces a value
     * each time it is called.
     */
    template <typename Generator>
    using generator_result_t = std::decay_t<std::invoke_result_t<Generator>>;

    /* Used to deduce the type of the numerator/denominator of a fraction. */
    template <typename Pair>
    struct pair_traits;

    template <typename T>
    struct pair_traits<std::pair<T, T>> {
        using value_type = T;
    };

    template <typename Pair>
    using pair_value_t = typename pair_traits<Pair>::value_type;

    /* Used to extract the "value type" of a complex type. */
    template <typename T>
    struct real_type {
        using type = T;
    };

    template <typename T>
    struct real_type<std::complex<T>> {
        using type = T;
    };

    template <typename T>
    using real_type_t = typename real_type<T>::type;

    // Return NaN, handling both real and complex types.
    template <typename T>
    SPECFUN_HOST_DEVICE inline std::enable_if_t<std::is_floating_point_v<T>, T> maybe_complex_NaN() {
        return std::numeric_limits<T>::quiet_NaN();
    }

    template <typename T>
    SPECFUN_HOST_DEVICE inline std::enable_if_t<!std::is_floating_point_v<T>, T> maybe_complex_NaN() {
        using V = typename T::value_type;
        return {std::numeric_limits<V>::quiet_NaN(), std::numeric_limits<V>::quiet_NaN()};
    }

    // Series evaluators.
    template <typename Generator, typename T = generator_result_t<Generator>>
    SPECFUN_HOST_DEVICE T series_eval(Generator &g, T init_val, real_type_t<T> tol, std::uint64_t max_terms,
                                      const char *func_name) {
        /* Sum an infinite series to a given precision.
         *
         * g : a generator of terms for the series.
         *
         * init_val : A starting value that terms are added to. This argument determines the
         *     type of the result.
         *
         * tol : relative tolerance for stopping criterion.
         *
         * max_terms : The maximum number of terms to add before giving up and declaring
         *     non-convergence.
         *
         * func_name : The name of the function within SciPy where this call to series_eval
         *     will ultimately be used. This is needed to pass to set_error in case
         *     of non-convergence.
         */
        T result = init_val;
        T term;
        for (std::uint64_t i = 0; i < max_terms; ++i) {
            term = g();
            result += term;
            if (std::abs(term) < std::abs(result) * tol) {
                return result;
            }
        }
        // Exceeded max terms without converging. Return NaN.
        set_error(func_name, SF_ERROR_NO_RESULT, NULL);
        return maybe_complex_NaN<T>();
    }

    template <typename Generator, typename T = generator_result_t<Generator>>
    SPECFUN_HOST_DEVICE T series_eval_fixed_length(Generator &g, T init_val, std::uint64_t num_terms) {
        /* Sum a fixed number of terms from a series.
         *
         * g : a generator of terms for the series.
         *
         * init_val : A starting value that terms are added to. This argument determines the
         *     type of the result.
         *
         * max_terms : The number of terms from the series to sum.
         *
         */
        T result = init_val;
        for (std::uint64_t i = 0; i < num_terms; ++i) {
            result += g();
        }
        return result;
    }

    /* Performs one step of Kahan summation. */
    template <typename T>
    SPECFUN_HOST_DEVICE void kahan_step(T& sum, T& comp, T x) {
        T y = x - comp;
        T t = sum + y;
        comp = (t - sum) - y;
        sum = t;
    }

    /* Evaluates an infinite series using Kahan summation.
     *
     * Denote the series by
     *
     *   S = a[0] + a[1] + a[2] + ...
     *
     * And for n = 0, 1, 2, ..., denote its n-th partial sum by
     *
     *   S[n] = a[0] + a[1] + ... + a[n]
     *
     * This function computes S[0], S[1], ... until a[n] is sufficiently
     * small or if the maximum number of terms have been evaluated.
     *
     * Parameters
     * ----------
     *   g
     *       Reference to generator that yields the sequence of values a[1],
     *       a[2], a[3], ...
     *
     *   tol
     *       Relative tolerance for convergence.  Specifically, stop iteration
     *       as soon as `abs(a[n]) <= tol * abs(S[n])` for some n >= 1.
     *
     *   max_terms
     *       Maximum number of terms after a[0] to evaluate.  It should be set
     *       large enough such that the convergence criterion is guaranteed
     *       to have been satisfied within that many terms if there is no
     *       rounding error.
     *
     *   init_val
     *       a[0].  Default is zero.  The type of this parameter (T) is used
     *       for intermediary computations as well as the result.
     *
     * Return Value
     * ------------
     * If the convergence criterion is satisfied by some `n <= max_terms`,
     * returns `(S[n], n)`.  Otherwise, returns `(S[max_terms], 0)`.
     */
    template <typename Generator, typename T = generator_result_t<Generator>>
    SPECFUN_HOST_DEVICE std::pair<T, std::uint64_t> series_eval_kahan(
        Generator &&g, real_type_t<T> tol, std::uint64_t max_terms, T init_val = T(0)) {

        T sum = init_val;
        T comp = 0;
        for (std::uint64_t i = 0; i < max_terms; ++i) {
            T term = g();
            kahan_step(sum, comp, term);
            if (std::abs(term) <= tol * std::abs(sum)) {
                return {sum, i + 1};
            }
        }
        return {sum, 0};
    }

    /* Generator that yields the difference of successive convergents of a
     * continued fraction.
     *
     * Let f[n] denote the n-th convergent of a continued fraction:
     *
     *                 a[1]   a[2]       a[n]
     *   f[n] = b[0] + ------ ------ ... ----
     *                 b[1] + b[2] +     b[n]
     *
     * with f[0] = b[0].  This generator yields the sequence of values
     * f[1]-f[0], f[2]-f[1], f[3]-f[2], ...
     *
     * Constructor Arguments
     * ---------------------
     *   cf
     *       Reference to generator that yields the terms of the continued
     *       fraction as (numerator, denominator) pairs, starting from
     *       (a[1], b[1]).
     *
     *       `cf` must outlive the ContinuedFractionSeriesGenerator object.
     *
     *       The constructed object always eagerly retrieves the next term
     *       of the continued fraction.  Specifically, (a[1], b[1]) is
     *       retrieved upon construction, and (a[n], b[n]) is retrieved after
     *       (n-1) calls of `()`.
     *
     * Type Arguments
     * --------------
     *   T
     *       Type in which computations are performed and results are turned.
     *
     * Remarks
     * -------
     * The series is computed using the recurrence relation described in [1].
     *
     * No error checking is performed.  The caller must ensure that all terms
     * are finite and that intermediary computations do not trigger floating
     * point exceptions such as overflow.
     *
     * The numerical stability of this method depends on the characteristics
     * of the continued fraction being evaluated.
     *
     * Reference
     * ---------
     * [1] Gautschi, W. (1967). “Computational Aspects of Three-Term
     *     Recurrence Relations.” SIAM Review, 9(1):24-82.
     */
    template <typename Generator, typename T = pair_value_t<generator_result_t<Generator>>>
    class ContinuedFractionSeriesGenerator {

    public:
        explicit ContinuedFractionSeriesGenerator(Generator &cf) : cf_(cf) {
            init();
        }

        double operator()() {
            double v = v_;
            advance();
            return v;
        }

    private:
        void init() {
            auto [num, denom] = cf_();
            T a = num;
            T b = denom;
            u_ = T(1);
            v_ = a / b;
            b_ = b;
        }

        void advance() {
            auto [num, denom] = cf_();
            T a = num;
            T b = denom;
            u_ = T(1) / (T(1) + (a * u_) / (b * b_));
            v_ *= (u_ - T(1));
            b_ = b;
        }

        Generator& cf_; // reference to continued fraction generator
        T v_;           // v[n] == f[n] - f[n-1], n >= 1
        T u_;           // u[1] = 1, u[n] = v[n]/v[n-1], n >= 2
        T b_;           // last denominator, i.e. b[n-1]
    };

    /* Converts a continued fraction into a series whose terms are the
     * difference of its successive convergents.
     *
     * See ContinuedFractionSeriesGenerator for details.
     */
    template <typename Generator, typename T = pair_value_t<generator_result_t<Generator>>>
    SPECFUN_HOST_DEVICE ContinuedFractionSeriesGenerator<Generator, T>
    continued_fraction_series(Generator &cf) {
        return ContinuedFractionSeriesGenerator<Generator, T>(cf);
    }

} // namespace detail
} // namespace special
