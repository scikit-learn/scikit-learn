#ifndef PYTHONIC_INCLUDE_NUMPY_NONZERO_HPP
#define PYTHONIC_INCLUDE_NUMPY_NONZERO_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E>
  auto nonzero(E const &expr)
      -> types::array_tuple<types::ndarray<long, types::array_tuple<long, 1>>, E::value>;

  DEFINE_FUNCTOR(pythonic::numpy, nonzero)
} // namespace numpy
PYTHONIC_NS_END

#endif
