#ifndef PYTHONIC_INCLUDE_NUMPY_DIFF_HPP
#define PYTHONIC_INCLUDE_NUMPY_DIFF_HPP

#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E>
  types::ndarray<typename E::dtype, types::array_tuple<long, E::value>>
  diff(E const &expr, long n = 1, long axis = -1);

  DEFINE_FUNCTOR(pythonic::numpy, diff);
} // namespace numpy
PYTHONIC_NS_END

#endif
