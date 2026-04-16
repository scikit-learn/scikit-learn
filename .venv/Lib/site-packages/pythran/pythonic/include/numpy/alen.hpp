#ifndef PYTHONIC_INCLUDE_NUMPY_ALEN_HPP
#define PYTHONIC_INCLUDE_NUMPY_ALEN_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  long alen(T &&expr);

  DEFINE_FUNCTOR(pythonic::numpy, alen);
} // namespace numpy
PYTHONIC_NS_END

#endif
