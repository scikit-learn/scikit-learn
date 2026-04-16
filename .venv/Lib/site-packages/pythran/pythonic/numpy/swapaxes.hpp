#ifndef PYTHONIC_NUMPY_SWAPAXES_HPP
#define PYTHONIC_NUMPY_SWAPAXES_HPP

#include "pythonic/include/numpy/swapaxes.hpp"

#include "pythonic/numpy/transpose.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  auto swapaxes(T &&a, int axis1, int axis2) -> decltype(functor::transpose{}(
      std::forward<T>(a), std::declval<types::array_tuple<long, std::decay_t<T>::value>>()))
  {
    constexpr long N = std::decay_t<T>::value;
    types::array_tuple<long, N> t;
    for (unsigned long i = 0; i < N; ++i)
      t[i] = i;
    std::swap(t[axis1], t[axis2]);
    return functor::transpose{}(std::forward<T>(a), t);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
