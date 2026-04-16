#ifndef PYTHONIC_NUMPY_RANDOM_DIRICHLET_HPP
#define PYTHONIC_NUMPY_RANDOM_DIRICHLET_HPP

#include "pythonic/include/numpy/random/dirichlet.hpp"
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
    types::ndarray<double, pS> dirichlet(double alpha, pS const &shape)
    {
      types::ndarray<double, pS> result{shape, types::none_type()};
      std::dirichlet_distribution<double> distribution{alpha};
      std::generate(result.fbegin(), result.fend(),
                    [&]() { return distribution(details::generator); });
      return result;
    }

    inline auto dirichlet(double alpha, long size)
        -> decltype(dirichlet(alpha, types::array_tuple<long, 1>{{size}}))
    {

      return dirichlet(alpha, types::array_tuple<long, 1>{{size}});
    }

    inline double dirichlet(double alpha, types::none_type d)
    {
      return std::dirichlet_distribution<double>{alpha}(details::generator);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
