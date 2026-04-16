#ifndef PYTHONIC_INCLUDE_NUMPY_ZEROSLIKE_HPP
#define PYTHONIC_INCLUDE_NUMPY_ZEROSLIKE_HPP

#include "pythonic/include/numpy/zeros.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E, class dtype>
  auto zeros_like(E const &expr, dtype d = dtype()) -> decltype(zeros(sutils::getshape(expr), d));

  template <class E>
  auto zeros_like(E const &expr, types::none_type d = builtins::None)
      -> decltype(zeros(sutils::getshape(expr), types::dtype_t<typename E::dtype>()));

  DEFINE_FUNCTOR(pythonic::numpy, zeros_like)
} // namespace numpy
PYTHONIC_NS_END

#endif
