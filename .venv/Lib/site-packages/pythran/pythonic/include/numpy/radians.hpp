#ifndef PYTHONIC_INCLUDE_NUMPY_RADIANS_HPP
#define PYTHONIC_INCLUDE_NUMPY_RADIANS_HPP

#include "pythonic/include/numpy/deg2rad.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  USING_FUNCTOR(radians, deg2rad);
}
PYTHONIC_NS_END

#endif
