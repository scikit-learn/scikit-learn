#ifndef PYTHONIC_NUMPY_RANDOM_RANDOM_HPP
#define PYTHONIC_NUMPY_RANDOM_RANDOM_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/random.hpp"

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
    types::ndarray<double, pS> random(pS const &shape)
    {
      types::ndarray<double, pS> result{shape, types::none_type()};
      std::uniform_real_distribution<double> distribution{0., 1.};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return distribution(details::generator); });
      return result;
    }

    inline auto random(long size) -> decltype(random(types::array_tuple<long, 1>{{size}}))
    {
      return random(types::array_tuple<long, 1>{{size}});
    }

    inline double random(types::none_type d)
    {
      return std::uniform_real_distribution<double>{0., 1.}(details::generator);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
