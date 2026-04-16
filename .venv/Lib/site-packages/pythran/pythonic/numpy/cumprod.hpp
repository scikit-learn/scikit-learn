#ifndef PYTHONIC_NUMPY_CUMPROD_HPP
#define PYTHONIC_NUMPY_CUMPROD_HPP

#include "pythonic/include/numpy/cumprod.hpp"

#include "pythonic/numpy/partial_sum.hpp"
#include "pythonic/operator_/imul.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E, class... Opts>
  auto cumprod(E &&e, Opts &&...opts)
      -> decltype(partial_sum<operator_::functor::imul>(std::forward<E>(e),
                                                        std::forward<Opts>(opts)...))
  {
    return partial_sum<operator_::functor::imul>(std::forward<E>(e), std::forward<Opts>(opts)...);
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(cumprod);
} // namespace numpy
PYTHONIC_NS_END

#endif
