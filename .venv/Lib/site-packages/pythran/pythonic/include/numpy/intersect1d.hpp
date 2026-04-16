#ifndef PYTHONIC_INCLUDE_NUMPY_INTERSECT1D_HPP
#define PYTHONIC_INCLUDE_NUMPY_INTERSECT1D_HPP

#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/types/combined.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class F>
  types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                 types::pshape<long>>
  intersect1d(E const &e, F const &f);

  DEFINE_FUNCTOR(pythonic::numpy, intersect1d);
} // namespace numpy
PYTHONIC_NS_END

#endif
