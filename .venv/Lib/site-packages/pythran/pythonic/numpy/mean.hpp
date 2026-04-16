#ifndef PYTHONIC_NUMPY_MEAN_HPP
#define PYTHONIC_NUMPY_MEAN_HPP

#include "pythonic/builtins/None.hpp"
#include "pythonic/include/numpy/mean.hpp"
#include "pythonic/numpy/asarray.hpp"
#include "pythonic/numpy/expand_dims.hpp"
#include "pythonic/numpy/sum.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E, class dtype>
  auto mean(E const &expr, types::none_type axis, dtype d, types::none_type out,
            types::false_immediate keepdims)
      -> decltype(sum(expr, axis, d) / details::dtype_or_double<dtype>(expr.flat_size()))
  {
    return sum(expr, axis, d) / details::dtype_or_double<dtype>(expr.flat_size());
  }

  template <class E, class dtype>
  auto mean(E const &expr, long axis, dtype d, types::none_type out,
            types::false_immediate keepdims) -> decltype(sum(expr, axis, d))
  {
    return sum(expr, axis, d) /= details::dtype_or_double<dtype>(sutils::getshape(expr)[axis]);
  }

  template <class E, class dtype>
  types::ndarray<details::dtype_or_double<dtype>,
                 typename details::make_scalar_pshape<E::value>::type>
  mean(E const &expr, types::none_type axis, dtype d, types::none_type out,
       types::true_immediate keep_dims)
  {
    return {typename details::make_scalar_pshape<E::value>::type(), mean(expr, axis, d, out)};
  }

  template <class E, class dtype>
  auto mean(E const &expr, long axis, dtype d, types::none_type out, types::true_immediate keepdims)
      -> decltype(expand_dims(mean(expr, axis, d), axis))
  {
    return expand_dims(mean(expr, axis, d), axis);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
