#ifndef PYTHONIC_INCLUDE_NUMPY_DEGREES_HPP
#define PYTHONIC_INCLUDE_NUMPY_DEGREES_HPP

#include "pythonic/include/numpy/rad2deg.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  USING_FUNCTOR(degrees, rad2deg);
}
PYTHONIC_NS_END

#endif
