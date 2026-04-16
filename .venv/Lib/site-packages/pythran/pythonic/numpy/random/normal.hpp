#ifndef PYTHONIC_NUMPY_RANDOM_NORMAL_HPP
#define PYTHONIC_NUMPY_RANDOM_NORMAL_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/normal.hpp"

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
    types::ndarray<double, pS> normal(double loc, double scale, pS const &shape)
    {
      types::ndarray<double, pS> result{shape, types::none_type()};
      std::normal_distribution<double> distribution{loc, scale};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return distribution(details::generator); });
      return result;
    }

    inline auto normal(double loc, double scale, long size)
        -> decltype(normal(loc, scale, types::array_tuple<long, 1>{{size}}))
    {
      return normal(loc, scale, types::array_tuple<long, 1>{{size}});
    }

    inline double normal(double loc, double scale, types::none_type d)
    {
      return std::normal_distribution<double>{loc, scale}(details::generator);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
