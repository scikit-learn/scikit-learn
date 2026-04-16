#ifndef PYTHONIC_NUMPY_MIN_HPP
#define PYTHONIC_NUMPY_MIN_HPP

#include "pythonic/include/numpy/min.hpp"

#include "pythonic/numpy/reduce.hpp"
#include "pythonic/operator_/imin.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class... Args>
  auto min(Args &&...args)
      -> decltype(reduce<operator_::functor::imin>(std::forward<Args>(args)...))
  {
    return reduce<operator_::functor::imin>(std::forward<Args>(args)...);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
