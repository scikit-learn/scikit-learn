#ifndef PYTHONIC_INCLUDE_NUMPY_REAL_HPP
#define PYTHONIC_INCLUDE_NUMPY_REAL_HPP

#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  auto real(E &&expr) -> decltype(builtins::getattr(types::attr::REAL{}, std::forward<E>(expr)));
  template <class T>
  auto real(types::list<T> const &expr) -> decltype(real(numpy::functor::asarray{}(expr)));

  DEFINE_FUNCTOR(pythonic::numpy, real);
} // namespace numpy
PYTHONIC_NS_END

#endif
