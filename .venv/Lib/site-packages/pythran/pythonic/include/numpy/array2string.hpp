#ifndef PYTHONIC_INCLUDE_NUMPY_ARRAY2STRING_HPP
#define PYTHONIC_INCLUDE_NUMPY_ARRAY2STRING_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::str array2string(E &&a);

  DEFINE_FUNCTOR(pythonic::numpy, array2string);
} // namespace numpy
PYTHONIC_NS_END

#endif
