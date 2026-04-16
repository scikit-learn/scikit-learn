#ifndef PYTHONIC_NUMPY_RANDOM_GAMMA_HPP
#define PYTHONIC_NUMPY_RANDOM_GAMMA_HPP

#include "pythonic/include/numpy/random/gamma.hpp"
#include "pythonic/include/numpy/random/generator.hpp"

#include "pythonic/types/NoneType.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"

#include <algorithm>
#include <random>

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {

    template <class pS>
    types::ndarray<double, pS> gamma(double shape, double scale, pS const &array_shape)
    {
      types::ndarray<double, pS> result{array_shape, types::none_type()};
      std::gamma_distribution<double> distribution{shape, scale};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return distribution(details::generator); });
      return result;
    }

    inline auto gamma(double shape, double scale, long size)
        -> decltype(gamma(shape, scale, types::array_tuple<long, 1>{{size}}))
    {
      return gamma(shape, scale, types::array_tuple<long, 1>{{size}});
    }

    inline double gamma(double shape, double scale, types::none_type d)
    {
      return std::gamma_distribution<double>{shape, scale}(details::generator);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
