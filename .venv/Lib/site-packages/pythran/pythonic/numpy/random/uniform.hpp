#ifndef PYTHONIC_NUMPY_RANDOM_RAYLEIGH_HPP
#define PYTHONIC_NUMPY_RANDOM_RAYLEIGH_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/uniform.hpp"

#include "pythonic/types/NoneType.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"
#include <math.h>

#include <algorithm>
#include <random>

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {

    template <class pS>
    types::ndarray<double, pS> uniform(double low, double high, pS const &array_shape)
    {
      types::ndarray<double, pS> result{array_shape, types::none_type()};
      std::generate(result.fbegin(), result.fend(), [&]() { return uniform(low, high); });
      return result;
    }

    inline auto uniform(double low, double high, long size)
        -> decltype(uniform(low, high, types::array_tuple<long, 1>{{size}}))
    {
      return uniform(low, high, types::array_tuple<long, 1>{{size}});
    }

    inline double uniform(double low, double high, types::none_type d)
    {
      return std::uniform_real_distribution<double>{low, high}(details::generator);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
