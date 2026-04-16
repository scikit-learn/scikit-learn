#ifndef PYTHONIC_INCLUDE_NUMPY_ARGMAX_HPP
#define PYTHONIC_INCLUDE_NUMPY_ARGMAX_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E>
  long argmax(E const &expr);

  template <class E>
  types::ndarray<long, types::array_tuple<long, E::value - 1>> argmax(E const &expr, long axis);

  DEFINE_FUNCTOR(pythonic::numpy, argmax);
} // namespace numpy
PYTHONIC_NS_END

#endif
