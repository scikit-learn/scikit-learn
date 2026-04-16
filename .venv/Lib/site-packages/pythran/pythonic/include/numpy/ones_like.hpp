#ifndef PYTHONIC_INCLUDE_NUMPY_ONESLIKE_HPP
#define PYTHONIC_INCLUDE_NUMPY_ONESLIKE_HPP

#include "pythonic/include/numpy/ones.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E, class dtype>
  auto ones_like(E const &expr, dtype d = dtype()) -> decltype(ones(sutils::getshape(expr), d));

  template <class E>
  auto ones_like(E const &expr, types::none_type d = builtins::None)
      -> decltype(ones(sutils::getshape(expr), types::dtype_t<typename E::dtype>()));

  DEFINE_FUNCTOR(pythonic::numpy, ones_like)
} // namespace numpy
PYTHONIC_NS_END

#endif
