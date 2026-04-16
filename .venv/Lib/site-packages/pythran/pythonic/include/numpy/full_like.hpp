#ifndef PYTHONIC_INCLUDE_NUMPY_FULLLIKE_HPP
#define PYTHONIC_INCLUDE_NUMPY_FULLLIKE_HPP

#include "pythonic/include/numpy/full.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E, class F, class dtype>
  auto full_like(E const &expr, F fill_value, dtype d = dtype())
      -> decltype(full(sutils::getshape(expr), fill_value, d));

  template <class E, class F>
  auto full_like(E const &expr, F fill_value, types::none_type d = builtins::None)
      -> decltype(full(sutils::getshape(expr), fill_value, types::dtype_t<typename E::dtype>()));

  DEFINE_FUNCTOR(pythonic::numpy, full_like)
} // namespace numpy
PYTHONIC_NS_END

#endif
