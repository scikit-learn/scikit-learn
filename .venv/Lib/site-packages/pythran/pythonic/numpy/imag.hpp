#ifndef PYTHONIC_NUMPY_IMAG_HPP
#define PYTHONIC_NUMPY_IMAG_HPP

#include "pythonic/include/numpy/imag.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/list.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  auto imag(E &&expr) -> decltype(builtins::getattr(types::attr::IMAG{}, std::forward<E>(expr)))
  {
    return builtins::getattr(types::attr::IMAG{}, std::forward<E>(expr));
  }

  template <class T>
  auto imag(types::list<T> const &expr) -> decltype(imag(numpy::functor::asarray{}(expr)))
  {
    return imag(numpy::functor::asarray{}(expr));
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
