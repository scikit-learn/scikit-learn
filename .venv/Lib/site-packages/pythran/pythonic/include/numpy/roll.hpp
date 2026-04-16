#ifndef PYTHONIC_INCLUDE_NUMPY_ROLL_HPP
#define PYTHONIC_INCLUDE_NUMPY_ROLL_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, pS> roll(types::ndarray<T, pS> const &expr, long shift);

  template <class T, class pS>
  types::ndarray<T, pS> roll(types::ndarray<T, pS> const &expr, long shift, long axis);

  template <class T, class pS, size_t N>
  types::ndarray<T, pS> roll(types::ndarray<T, pS> const &expr, types::array_tuple<long, N> shift,
                             types::array_tuple<long, N> axis);

  NUMPY_EXPR_TO_NDARRAY0_DECL(roll);
  DEFINE_FUNCTOR(pythonic::numpy, roll);
} // namespace numpy
PYTHONIC_NS_END

#endif
