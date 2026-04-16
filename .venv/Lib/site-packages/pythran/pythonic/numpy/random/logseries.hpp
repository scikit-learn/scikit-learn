#ifndef PYTHONIC_NUMPY_RANDOM_LOGSERIES_HPP
#define PYTHONIC_NUMPY_RANDOM_LOGSERIES_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/logseries.hpp"

#include "pythonic/types/NoneType.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"
#include <math.h>

#include <algorithm>
#include <random>

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {

    template <class pS>
    types::ndarray<double, pS> logseries(double p, pS const &shape)
    {
      types::ndarray<double, pS> result{shape, types::none_type()};
      std::generate(result.fbegin(), result.fend(), [&]() { return logseries(p); });
      return result;
    }

    inline auto logseries(double p, long size)
        -> decltype(logseries(p, types::array_tuple<long, 1>{{size}}))
    {
      return logseries(p, types::array_tuple<long, 1>{{size}});
    }

    inline double logseries(double p, types::none_type d)
    {
      double q, r, U, V;
      double result;

      r = log1p(-p);

      while (1) {
        V = std::uniform_real_distribution<double>{0., 1.}(details::generator);
        if (V >= p) {
          return 1;
        }
        U = std::uniform_real_distribution<double>{0., 1.}(details::generator);
        q = -expm1(r * U);
        if (V <= q * q) {
          result = (double)floor(1 + log(V) / log(q));
          if ((result < 1) || (V == 0.0)) {
            continue;
          } else {
            return result;
          }
        }
        if (V >= q) {
          return 1;
        }
        return 2;
      }
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
