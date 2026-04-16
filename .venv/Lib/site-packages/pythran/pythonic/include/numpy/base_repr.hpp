#ifndef PYTHONIC_INCLUDE_NUMPY_BASEREPR_HPP
#define PYTHONIC_INCLUDE_NUMPY_BASEREPR_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  types::str base_repr(long number, long base = 2, long padding = 0);

  DEFINE_FUNCTOR(pythonic::numpy, base_repr);
} // namespace numpy
PYTHONIC_NS_END

#endif
