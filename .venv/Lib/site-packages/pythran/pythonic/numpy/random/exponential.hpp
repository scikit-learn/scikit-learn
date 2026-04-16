#ifndef PYTHONIC_NUMPY_RANDOM_EXPONENTIAL_HPP
#define PYTHONIC_NUMPY_RANDOM_EXPONENTIAL_HPP

#include "pythonic/include/numpy/random/exponential.hpp"
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
    types::ndarray<double, pS> exponential(double scale, pS const &shape)
    {
      types::ndarray<double, pS> result{shape, types::none_type()};
      std::exponential_distribution<double> distribution{1 / scale};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return distribution(details::generator); });
      return result;
    }

    inline auto exponential(double scale, long size)
        -> decltype(exponential(scale, types::array_tuple<long, 1>{{size}}))
    {

      return exponential(scale, types::array_tuple<long, 1>{{size}});
    }

    inline double exponential(double scale, types::none_type d)
    {
      return std::exponential_distribution<double>{1 / scale}(details::generator);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
