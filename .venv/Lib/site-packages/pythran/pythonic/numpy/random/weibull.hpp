#ifndef PYTHONIC_NUMPY_RANDOM_WEIBULL_HPP
#define PYTHONIC_NUMPY_RANDOM_WEIBULL_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/weibull.hpp"

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
    types::ndarray<double, pS> weibull(double a, pS const &shape)
    {
      types::ndarray<double, pS> result{shape, types::none_type()};
      std::weibull_distribution<double> distribution{a};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return distribution(details::generator); });
      return result;
    }

    inline auto weibull(double a, long size)
        -> decltype(weibull(a, types::array_tuple<long, 1>{{size}}))
    {

      return weibull(a, types::array_tuple<long, 1>{{size}});
    }

    inline double weibull(double a, types::none_type d)
    {
      return std::weibull_distribution<double>{a}(details::generator);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
