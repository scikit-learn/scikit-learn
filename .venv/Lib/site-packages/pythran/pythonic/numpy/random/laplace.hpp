#ifndef PYTHONIC_NUMPY_RANDOM_LAPLACE_HPP
#define PYTHONIC_NUMPY_RANDOM_LAPLACE_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/laplace.hpp"

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
    types::ndarray<double, pS> laplace(double loc, double scale, pS const &shape)
    {
      types::ndarray<double, pS> result{shape, types::none_type()};
      std::generate(result.fbegin(), result.fend(), [&]() { return laplace(loc, scale); });
      return result;
    }

    inline auto laplace(double loc, double scale, long size)
        -> decltype(laplace(loc, scale, types::array_tuple<long, 1>{{size}}))
    {
      return laplace(loc, scale, types::array_tuple<long, 1>{{size}});
    }

    inline double laplace(double loc, double scale, types::none_type d)
    {
      double U = std::uniform_real_distribution<double>{0., 1.}(details::generator);
      if (U >= 0.5) {
        U = loc - scale * xsimd::log(2.0 - U - U);
      } else if (U > 0.0) {
        U = loc + scale * xsimd::log(U + U);
      } else {
        U = laplace(loc, scale);
      }
      return U;
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
