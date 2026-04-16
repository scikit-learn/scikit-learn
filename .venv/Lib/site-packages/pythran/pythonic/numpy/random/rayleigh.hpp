#ifndef PYTHONIC_NUMPY_RANDOM_RAYLEIGH_HPP
#define PYTHONIC_NUMPY_RANDOM_RAYLEIGH_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/rayleigh.hpp"

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
    types::ndarray<double, pS> rayleigh(double scale, pS const &array_shape)
    {
      types::ndarray<double, pS> result{array_shape, types::none_type()};
      std::generate(result.fbegin(), result.fend(), [&]() { return rayleigh(scale); });
      return result;
    }

    inline auto rayleigh(double scale, long size)
        -> decltype(rayleigh(scale, types::array_tuple<long, 1>{{size}}))
    {
      return rayleigh(scale, types::array_tuple<long, 1>{{size}});
    }

    inline double rayleigh(double scale, types::none_type d)
    {
      return scale * sqrt(-2.0 * log(1.0 - std::uniform_real_distribution<double>{0., 1.}(
                                               details::generator)));
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
