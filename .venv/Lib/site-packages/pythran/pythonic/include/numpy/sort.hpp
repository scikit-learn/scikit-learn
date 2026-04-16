#ifndef PYTHONIC_INCLUDE_NUMPY_SORT_HPP
#define PYTHONIC_INCLUDE_NUMPY_SORT_HPP

#include <algorithm>

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  template <class E>
  types::ndarray<typename E::dtype, types::array_tuple<long, 1>> sort(E const &expr,
                                                                      types::none_type);

  template <class E>
  types::ndarray<typename E::dtype, types::array_tuple<long, E::value>> sort(E const &expr,
                                                                             long axis = -1);

  template <class E>
  types::ndarray<typename E::dtype, types::array_tuple<long, E::value>>
  sort(E const &expr, long axis, types::str const &kind);

  NUMPY_EXPR_TO_NDARRAY0_DECL(sort);
  DEFINE_FUNCTOR(pythonic::numpy, sort);
} // namespace numpy
PYTHONIC_NS_END

#endif
