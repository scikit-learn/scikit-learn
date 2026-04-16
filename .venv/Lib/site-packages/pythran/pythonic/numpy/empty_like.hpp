#ifndef PYTHONIC_NUMPY_EMPTYLIKE_HPP
#define PYTHONIC_NUMPY_EMPTYLIKE_HPP

#include "pythonic/include/numpy/empty_like.hpp"

#include "pythonic/numpy/empty.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class dtype>
  auto empty_like(E const &expr, dtype d) -> decltype(empty(sutils::getshape(expr), d))
  {
    return empty(sutils::getshape(expr), d);
  }

  template <class E>
  auto empty_like(E const &expr, types::none_type)
      -> decltype(empty(sutils::getshape(expr), types::dtype_t<typename E::dtype>()))
  {
    return empty(sutils::getshape(expr), types::dtype_t<typename E::dtype>());
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
