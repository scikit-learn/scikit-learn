#ifndef PYTHONIC_INCLUDE_NUMPY_DIAG_HPP
#define PYTHONIC_INCLUDE_NUMPY_DIAG_HPP

#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  std::enable_if_t<std::tuple_size<pS>::value == 2, types::ndarray<T, types::pshape<long>>>
  diag(types::ndarray<T, pS> const &a, long k = 0);

  template <class T, class pS>
  std::enable_if_t<std::tuple_size<pS>::value == 1, types::ndarray<T, types::array_tuple<long, 2>>>
  diag(types::ndarray<T, pS> const &a, long k = 0);

  template <class T>
  auto diag(types::list<T> const &a, long k = 0) -> decltype(diag(asarray(a), k));

  NUMPY_EXPR_TO_NDARRAY0_DECL(diag);
  DEFINE_FUNCTOR(pythonic::numpy, diag);
} // namespace numpy
PYTHONIC_NS_END

#endif
