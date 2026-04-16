#ifndef PYTHONIC_NUMPY_ZEROSLIKE_HPP
#define PYTHONIC_NUMPY_ZEROSLIKE_HPP

#include "pythonic/include/numpy/zeros_like.hpp"

#include "pythonic/numpy/zeros.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E, class dtype>
  auto zeros_like(E const &expr, dtype d) -> decltype(zeros(sutils::getshape(expr), d))
  {
    return zeros(sutils::getshape(expr), d);
  }

  template <class E>
  auto zeros_like(E const &expr, types::none_type)
      -> decltype(zeros(sutils::getshape(expr), types::dtype_t<typename E::dtype>()))
  {
    return zeros(sutils::getshape(expr), types::dtype_t<typename E::dtype>());
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
