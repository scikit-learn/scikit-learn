#ifndef PYTHONIC_INCLUDE_NUMPY_IMAG_HPP
#define PYTHONIC_INCLUDE_NUMPY_IMAG_HPP

#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  auto imag(E &&expr) -> decltype(builtins::getattr(types::attr::IMAG{}, std::forward<E>(expr)));

  template <class T>
  auto imag(types::list<T> const &expr) -> decltype(imag(numpy::functor::asarray{}(expr)));

  DEFINE_FUNCTOR(pythonic::numpy, imag);
} // namespace numpy
PYTHONIC_NS_END

#endif
