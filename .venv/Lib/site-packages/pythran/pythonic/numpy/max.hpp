#ifndef PYTHONIC_NUMPY_MAX_HPP
#define PYTHONIC_NUMPY_MAX_HPP

#include "pythonic/include/numpy/max.hpp"

#include "pythonic/numpy/reduce.hpp"
#include "pythonic/operator_/imax.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class... Args>
  auto max(Args &&...args)
      -> decltype(reduce<operator_::functor::imax>(std::forward<Args>(args)...))
  {
    return reduce<operator_::functor::imax>(std::forward<Args>(args)...);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
