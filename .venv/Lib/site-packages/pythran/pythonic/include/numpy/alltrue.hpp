#ifndef PYTHONIC_INCLUDE_NUMPY_ALLTRUE_HPP
#define PYTHONIC_INCLUDE_NUMPY_ALLTRUE_HPP

#include "pythonic/include/numpy/all.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class... Types>
  auto alltrue(Types &&...types) -> decltype(all(std::forward<Types>(types)...));

  DEFINE_FUNCTOR(pythonic::numpy, alltrue);
} // namespace numpy
PYTHONIC_NS_END

#endif
