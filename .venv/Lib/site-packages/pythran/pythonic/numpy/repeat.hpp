#ifndef PYTHONIC_NUMPY_REPEAT_HPP
#define PYTHONIC_NUMPY_REPEAT_HPP

#include "pythonic/include/numpy/repeat.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
  repeat(types::ndarray<T, pS> const &expr, long repeats, long axis)
  {
    constexpr auto N = std::tuple_size<pS>::value;
    if (axis < 0)
      axis += N;

    auto shape = sutils::getshape(expr);
    const long stride =
        std::accumulate(shape.begin() + axis + 1, shape.end(), 1L, std::multiplies<long>());
    shape[axis] *= repeats;

    types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>> out(shape,
                                                                                builtins::None);
    auto out_iter = out.fbegin();
    for (auto iter = expr.fbegin(), end = expr.fend(); iter != end; iter += stride)
      for (int i = 0; i < repeats; ++i)
        out_iter = std::copy(iter, iter + stride, out_iter);
    return out;
  }
  template <class T, class pS>
  types::ndarray<T, types::pshape<long>> repeat(types::ndarray<T, pS> const &expr, long repeats,
                                                types::none_type axis)
  {
    types::ndarray<T, types::pshape<long>> out(types::pshape<long>{expr.flat_size() * repeats},
                                               builtins::None);
    auto out_iter = out.fbegin();
    for (auto iter = expr.fbegin(), end = expr.fend(); iter != end; ++iter)
      for (int i = 0; i < repeats; ++i)
        *out_iter++ = *iter;
    return out;
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(repeat);
} // namespace numpy
PYTHONIC_NS_END

#endif
