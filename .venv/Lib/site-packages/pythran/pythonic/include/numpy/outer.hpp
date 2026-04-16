#ifndef PYTHONIC_INCLUDE_NUMPY_OUTER_HPP
#define PYTHONIC_INCLUDE_NUMPY_OUTER_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T0, class pS0, class T1, class pS1>
  types::ndarray<decltype(std::declval<T0>() + std::declval<T1>()), types::pshape<long, long>>
  outer(types::ndarray<T0, pS0> const &a, types::ndarray<T1, pS1> const &b);

  template <class T0, class pS0, class E1>
  auto outer(types::ndarray<T0, pS0> const &a, E1 const &b) -> decltype(outer(a, asarray(b)));

  template <class E0, class T1, class pS1>
  auto outer(E0 const &a, types::ndarray<T1, pS1> const &b) -> decltype(outer(asarray(a), b));

  template <class E0, class E1>
  auto outer(E0 const &a, E1 const &b) -> decltype(outer(asarray(a), asarray(b)));

  DEFINE_FUNCTOR(pythonic::numpy, outer);
} // namespace numpy
PYTHONIC_NS_END

#endif
