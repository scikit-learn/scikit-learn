#ifndef PYTHONIC_NUMPY_RANDOM_PARETO_HPP
#define PYTHONIC_NUMPY_RANDOM_PARETO_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/pareto.hpp"

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
    types::ndarray<double, pS> pareto(double a, pS const &shape)
    {
      types::ndarray<double, pS> result{shape, types::none_type()};
      std::exponential_distribution<double> distribution{};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return expm1(distribution(details::generator) / a); });
      return result;
    }

    inline auto pareto(double a, long size)
        -> decltype(pareto(a, types::array_tuple<long, 1>{{size}}))
    {

      return pareto(a, types::array_tuple<long, 1>{{size}});
    }

    inline double pareto(double a, types::none_type d)
    {
      return expm1(std::exponential_distribution<double>{}(details::generator) / a);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
