#ifndef PYTHONIC_INCLUDE_NUMPY_CONJ_HPP
#define PYTHONIC_INCLUDE_NUMPY_CONJ_HPP

#include "pythonic/include/numpy/conjugate.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  USING_FUNCTOR(conj, numpy::functor::conjugate);
}
PYTHONIC_NS_END

#endif
