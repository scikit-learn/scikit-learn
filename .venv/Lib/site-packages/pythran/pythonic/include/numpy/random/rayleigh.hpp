#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_RAYLEIGH_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_RAYLEIGH_HPP

#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"
#include <math.h>

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    template <class pS>
    types::ndarray<double, pS> rayleigh(double scale, pS const &array_shape);

    auto rayleigh(double scale, long size)
        -> decltype(rayleigh(scale, types::array_tuple<long, 1>{{size}}));

    double rayleigh(double scale = 1.0, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, rayleigh);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
