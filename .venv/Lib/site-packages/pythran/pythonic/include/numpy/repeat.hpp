#ifndef PYTHONIC_INCLUDE_NUMPY_REPEAT_HPP
#define PYTHONIC_INCLUDE_NUMPY_REPEAT_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
  repeat(types::ndarray<T, pS> const &expr, long repeats, long axis);

  template <class T, class pS>
  types::ndarray<T, types::pshape<long>> repeat(types::ndarray<T, pS> const &expr, long repeats,
                                                types::none_type axis = types::none_type{});

  NUMPY_EXPR_TO_NDARRAY0_DECL(repeat);
  DEFINE_FUNCTOR(pythonic::numpy, repeat);
} // namespace numpy
PYTHONIC_NS_END

#endif
