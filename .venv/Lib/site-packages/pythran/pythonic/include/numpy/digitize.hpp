#ifndef PYTHONIC_INCLUDE_NUMPY_DIGITIZE_HPP
#define PYTHONIC_INCLUDE_NUMPY_DIGITIZE_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/operator_/gt.hpp"
#include "pythonic/include/operator_/lt.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class F>
  types::ndarray<long, types::pshape<long>> digitize(E const &expr, F const &b);

  DEFINE_FUNCTOR(pythonic::numpy, digitize);
} // namespace numpy
PYTHONIC_NS_END

#endif
