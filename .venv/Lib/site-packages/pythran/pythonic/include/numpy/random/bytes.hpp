#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_BYTES_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_BYTES_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    types::str bytes(long length);

    DEFINE_FUNCTOR(pythonic::numpy::random, bytes);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
