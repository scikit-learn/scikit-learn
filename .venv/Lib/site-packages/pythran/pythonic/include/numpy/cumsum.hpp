#ifndef PYTHONIC_INCLUDE_NUMPY_CUMSUM_HPP
#define PYTHONIC_INCLUDE_NUMPY_CUMSUM_HPP

#include "pythonic/include/numpy/partial_sum.hpp"
#include "pythonic/include/operator_/iadd.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E, class... Opts>
  auto cumsum(E &&e, Opts &&...opts)
      -> decltype(partial_sum<operator_::functor::add>(std::forward<E>(e),
                                                       std::forward<Opts>(opts)...));

  DEFINE_FUNCTOR(pythonic::numpy, cumsum);
} // namespace numpy
PYTHONIC_NS_END

#endif
