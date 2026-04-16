#ifndef PYTHONIC_NUMPY_ALLTRUE_HPP
#define PYTHONIC_NUMPY_ALLTRUE_HPP

#include "pythonic/include/numpy/alltrue.hpp"

#include "pythonic/numpy/all.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class... Types>
  auto alltrue(Types &&...types) -> decltype(all(std::forward<Types>(types)...))
  {
    return all(std::forward<Types>(types)...);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
