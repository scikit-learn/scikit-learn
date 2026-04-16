#ifndef PYTHONIC_INCLUDE_NUMPY_SUM_HPP
#define PYTHONIC_INCLUDE_NUMPY_SUM_HPP

#include "pythonic/include/numpy/reduce.hpp"
#include "pythonic/include/operator_/iadd.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class... Args>
  auto sum(Args &&...args)
      -> decltype(reduce<operator_::functor::iadd>(std::forward<Args>(args)...))
  {
    return reduce<operator_::functor::iadd>(std::forward<Args>(args)...);
  }

  DEFINE_FUNCTOR(pythonic::numpy, sum);
} // namespace numpy
PYTHONIC_NS_END

#endif
