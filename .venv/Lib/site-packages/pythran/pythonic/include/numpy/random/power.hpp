#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_POWER_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_POWER_HPP

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
    types::ndarray<double, pS> power(double a, pS const &shape);
    auto power(double a, long size) -> decltype(power(a, types::array_tuple<long, 1>{{size}}));

    double power(double a, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, power);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
