#ifndef PYTHONIC_INCLUDE_NUMPY_FABS_HPP
#define PYTHONIC_INCLUDE_NUMPY_FABS_HPP

#include "pythonic/include/numpy/abs.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  USING_FUNCTOR(fabs, numpy::functor::abs);
}
PYTHONIC_NS_END

#endif
