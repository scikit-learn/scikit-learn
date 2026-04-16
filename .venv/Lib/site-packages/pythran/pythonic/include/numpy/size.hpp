#ifndef PYTHONIC_INCLUDE_NUMPY_SIZE_HPP
#define PYTHONIC_INCLUDE_NUMPY_SIZE_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E>
  auto size(E const &e) -> decltype(e.flat_size());

  inline long size(...)
  {
    return 1;
  }

  DEFINE_FUNCTOR(pythonic::numpy, size)
} // namespace numpy
PYTHONIC_NS_END

#endif
