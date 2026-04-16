#ifndef PYTHONIC_INCLUDE_NUMPY_MAX_HPP
#define PYTHONIC_INCLUDE_NUMPY_MAX_HPP

#include "pythonic/include/numpy/reduce.hpp"
#include "pythonic/include/operator_/imax.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class... Args>
  auto max(Args &&...args)
      -> decltype(reduce<operator_::functor::imax>(std::forward<Args>(args)...));

  DEFINE_FUNCTOR(pythonic::numpy, max);
} // namespace numpy
PYTHONIC_NS_END

#endif
