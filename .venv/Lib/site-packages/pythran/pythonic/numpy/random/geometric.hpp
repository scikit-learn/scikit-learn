#ifndef PYTHONIC_NUMPY_RANDOM_GEOMETRIC_HPP
#define PYTHONIC_NUMPY_RANDOM_GEOMETRIC_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/geometric.hpp"

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
    types::ndarray<double, pS> geometric(double p, pS const &shape)
    {
      types::ndarray<double, pS> result{shape, types::none_type()};
      std::geometric_distribution<int> distribution{p};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return distribution(details::generator); });
      return result;
    }

    inline auto geometric(double p, long size)
        -> decltype(geometric(p, types::array_tuple<long, 1>{{size}}))
    {
      return geometric(p, types::array_tuple<long, 1>{{size}});
    }

    inline double geometric(double p, types::none_type d)
    {
      return std::geometric_distribution<int>{p}(details::generator);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
