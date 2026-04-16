#ifndef PYTHONIC_SCIPY_SPECIAL_BINOM_HPP
#define PYTHONIC_SCIPY_SPECIAL_BINOM_HPP

#include "pythonic/include/scipy/special/binom.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

#include "pythonic/utils/boost_local_config.hpp"
#include <boost/math/special_functions/binomial.hpp>

PYTHONIC_NS_BEGIN

namespace scipy
{
  namespace special
  {
    namespace details
    {
      template <class T0, class T1>
      double binom(T0 n, T1 k)
      {
        static_assert(std::is_integral<T0>::value && std::is_integral<T1>::value,
                      "only support integer case of scipy.special.binom");
        using namespace boost::math::policies;
        return boost::math::binomial_coefficient<double>(n, k, make_policy(promote_double<true>()));
      }
    } // namespace details

#define NUMPY_NARY_FUNC_NAME binom
#define NUMPY_NARY_FUNC_SYM details::binom
#include "pythonic/types/numpy_nary_expr.hpp"
  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
