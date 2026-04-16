#ifndef PYTHONIC_INCLUDE_NUMPY_ARGSORT_HPP
#define PYTHONIC_INCLUDE_NUMPY_ARGSORT_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::ndarray<long, types::array_tuple<long, 1>> argsort(E const &expr, types::none_type,
                                                            types::none_type = {});

  template <class T, class pS>
  types::ndarray<long, pS> argsort(types::ndarray<T, pS> const &a, long axis = -1,
                                   types::none_type kind = {});

  template <class T, class pS>
  types::ndarray<long, pS> argsort(types::ndarray<T, pS> const &a, long axis,
                                   types::str const &kind);

  NUMPY_EXPR_TO_NDARRAY0_DECL(argsort);

  DEFINE_FUNCTOR(pythonic::numpy, argsort);
} // namespace numpy
PYTHONIC_NS_END

#endif
