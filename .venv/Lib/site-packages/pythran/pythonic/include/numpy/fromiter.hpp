#ifndef PYTHONIC_INCLUDE_NUMPY_FROMITER_HPP
#define PYTHONIC_INCLUDE_NUMPY_FROMITER_HPP

#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class Iterable, class dtype = functor::float64>
  types::ndarray<typename std::remove_cv_t<std::remove_reference_t<Iterable>>::value_type,
                 types::pshape<long>>
  fromiter(Iterable &&iterable, dtype d = dtype(), long count = -1);

  DEFINE_FUNCTOR(pythonic::numpy, fromiter);
} // namespace numpy
PYTHONIC_NS_END

#endif
