#ifndef PYTHONIC_NUMPY_RANDOM_NEGATIVE_BINOMIAL_HPP
#define PYTHONIC_NUMPY_RANDOM_NEGATIVE_BINOMIAL_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/negative_binomial.hpp"

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
    types::ndarray<long, pS> negative_binomial(long n, double p, pS const &shape)
    {
      types::ndarray<long, pS> result{shape, types::none_type()};
      std::negative_binomial_distribution<long> distribution{n, p};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return distribution(details::generator); });
      return result;
    }

    inline auto negative_binomial(long n, double p, long size)
        -> decltype(negative_binomial(n, p, types::array_tuple<long, 1>{{size}}))
    {
      return negative_binomial(n, p, types::array_tuple<long, 1>{{size}});
    }

    inline long negative_binomial(long n, double p, types::none_type d)
    {
      std::negative_binomial_distribution<long> distribution{n, p};
      return distribution(details::generator);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
