#ifndef PYTHONIC_NUMPY_SIZE_HPP
#define PYTHONIC_NUMPY_SIZE_HPP

#include "pythonic/include/numpy/size.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E>
  auto size(E const &e) -> decltype(e.flat_size())
  {
    return e.flat_size();
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
