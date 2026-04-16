#ifndef PYTHONIC_NUMPY_RANDOM_LOGNORMAL_HPP
#define PYTHONIC_NUMPY_RANDOM_LOGNORMAL_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/lognormal.hpp"

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
    types::ndarray<double, pS> lognormal(double mean, double sigma, pS const &shape)
    {
      types::ndarray<double, pS> result{shape, types::none_type()};
      std::lognormal_distribution<double> distribution{mean, sigma};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return distribution(details::generator); });
      return result;
    }

    inline auto lognormal(double mean, double sigma, long size)
        -> decltype(lognormal(mean, sigma, types::array_tuple<long, 1>{{size}}))
    {
      return lognormal(mean, sigma, types::array_tuple<long, 1>{{size}});
    }

    inline double lognormal(double mean, double sigma, types::none_type d)
    {
      return std::lognormal_distribution<double>{mean, sigma}(details::generator);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
