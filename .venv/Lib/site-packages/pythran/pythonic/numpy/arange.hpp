#ifndef PYTHONIC_NUMPY_ARANGE_HPP
#define PYTHONIC_NUMPY_ARANGE_HPP

#include "pythonic/include/numpy/arange.hpp"

#include "pythonic/operator_/pos.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class U, class S, class dtype>
  types::numpy_expr<pythonic::operator_::functor::pos, details::arange_index<typename dtype::type>>
  arange(T begin, U end, S step, dtype d)
  {
    using R = typename dtype::type;
    long size;
    if (std::is_integral<R>::value)
      size = std::max(R(0), R((end - begin + step - 1) / step));
    else
      size = std::max(R(0), R(std::ceil((end - begin) / step)));
    return {details::arange_index<R>{(R)begin, (R)step, size}};
  }

  template <class T>
  types::numpy_expr<pythonic::operator_::functor::pos,
                    details::arange_index<typename types::dtype_t<T>::type>>
  arange(T end)
  {
    return arange<T, T, T, types::dtype_t<T>>(T(0), end);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
