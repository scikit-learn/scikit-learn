#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_UNIFORM_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_UNIFORM_HPP

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
    types::ndarray<double, pS> uniform(double low, double high, pS const &array_shape);

    auto uniform(double low, double high, long size)
        -> decltype(uniform(low, high, types::array_tuple<long, 1>{{size}}));

    double uniform(double low = 0.0, double high = 1.0, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, uniform);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
