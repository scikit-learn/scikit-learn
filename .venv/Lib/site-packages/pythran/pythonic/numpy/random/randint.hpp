#ifndef PYTHONIC_NUMPY_RANDOM_RANDINT_HPP
#define PYTHONIC_NUMPY_RANDOM_RANDINT_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/randint.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"

#include <random>

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {

    template <class pS>
    std::enable_if_t<!std::is_integral<pS>::value, types::ndarray<long, pS>>
    randint(long min, long max, pS const &shape)
    {
      types::ndarray<long, pS> result{shape, types::none_type()};
      std::uniform_int_distribution<long> distribution{min, max - 1};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return distribution(details::generator); });
      return result;
    }

    template <class pS>
    std::enable_if_t<std::is_integral<pS>::value, types::ndarray<long, types::pshape<long>>>
    randint(long min, long max, pS const &shape)
    {
      return randint(min, max, types::pshape<long>{shape});
    }

    template <class pS>
    auto randint(long max, types::none_type, pS const &shape) -> decltype(randint(0, max, shape))
    {
      return randint(0, max, shape);
    }

    inline auto randint(long min, long max, long size)
        -> decltype(randint(min, max, types::array_tuple<long, 1>{{size}}))
    {
      return randint(min, max, types::array_tuple<long, 1>{{size}});
    }

    inline long randint(long max, types::none_type)
    {
      return std::uniform_int_distribution<long>{0, max - 1}(details::generator);
    }

    inline long randint(long min, long max)
    {
      return std::uniform_int_distribution<long>{min, max - 1}(details::generator);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
