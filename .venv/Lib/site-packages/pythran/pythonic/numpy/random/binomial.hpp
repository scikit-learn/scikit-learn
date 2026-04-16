#ifndef PYTHONIC_NUMPY_RANDOM_BINOMIAL_HPP
#define PYTHONIC_NUMPY_RANDOM_BINOMIAL_HPP

#include "pythonic/include/numpy/random/binomial.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/numpy_expr.hpp"

#include "pythonic/types/exceptions.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    namespace details
    {
      inline void parameters_check(double n, double p)
      {
        if (n < 0)
          throw pythonic::types::ValueError("n < 0");
        if (p < 0)
          throw pythonic::types::ValueError("p < 0");
        else if (p > 1)
          throw pythonic::types::ValueError("p > 1");
      }
    } // namespace details

    template <class pS>
    types::ndarray<long, pS> binomial(double n, double p, pS const &shape)
    {
      details::parameters_check(n, p);
      types::ndarray<long, pS> result{shape, types::none_type()};
      std::binomial_distribution<long> distribution{(long)n, p};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return distribution(details::generator); });
      return result;
    }

    inline auto binomial(double n, double p, long size)
        -> decltype(binomial(n, p, types::array_tuple<long, 1>{{size}}))
    {
      return binomial(n, p, types::array_tuple<long, 1>{{size}});
    }

    inline long binomial(double n, double p, types::none_type d)
    {
      details::parameters_check(n, p);
      return std::binomial_distribution<long>{(long)n, p}(details::generator);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
