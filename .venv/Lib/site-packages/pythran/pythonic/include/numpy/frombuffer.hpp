#ifndef PYTHONIC_INCLUDE_NUMPY_FROMBUFFER_HPP
#define PYTHONIC_INCLUDE_NUMPY_FROMBUFFER_HPP

#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <limits>
#include <sstream>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::pshape<long>>
  frombuffer(types::str const &string, dtype d = dtype(), long count = -1, long offset = 0);

  DEFINE_FUNCTOR(pythonic::numpy, frombuffer);
} // namespace numpy
PYTHONIC_NS_END

#endif
