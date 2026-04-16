#ifndef PYTHONIC_INCLUDE_NUMPY_EMPTYLIKE_HPP
#define PYTHONIC_INCLUDE_NUMPY_EMPTYLIKE_HPP

#include "pythonic/include/numpy/empty.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class dtype>
  auto empty_like(E const &expr, dtype d = dtype()) -> decltype(empty(sutils::getshape(expr), d));

  template <class E>
  auto empty_like(E const &expr, types::none_type d = builtins::None)
      -> decltype(empty(sutils::getshape(expr), types::dtype_t<typename E::dtype>()));

  DEFINE_FUNCTOR(pythonic::numpy, empty_like)
} // namespace numpy
PYTHONIC_NS_END

#endif
