#ifndef PYTHONIC_INCLUDE_NUMPY_TILE_HPP
#define PYTHONIC_INCLUDE_NUMPY_TILE_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::ndarray<typename E::dtype, types::array_tuple<long, E::value>> tile(E const &expr,
                                                                             long reps);

  template <class E, size_t N>
  types::ndarray<typename E::dtype, types::array_tuple<long, N>>
  tile(E const &expr, types::array_tuple<long, N> const &reps);

  DEFINE_FUNCTOR(pythonic::numpy, tile);
} // namespace numpy
PYTHONIC_NS_END

#endif
