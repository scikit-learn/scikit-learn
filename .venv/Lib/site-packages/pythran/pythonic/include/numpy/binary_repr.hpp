#ifndef PYTHONIC_INCLUDE_NUMPY_BINARYREPR_HPP
#define PYTHONIC_INCLUDE_NUMPY_BINARYREPR_HPP

#include "pythonic/include/numpy/base_repr.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  types::str binary_repr(long number, types::none_type width = builtins::None);

  types::str binary_repr(long number, long width);

  DEFINE_FUNCTOR(pythonic::numpy, binary_repr);
} // namespace numpy
PYTHONIC_NS_END

#endif
