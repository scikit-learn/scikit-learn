#ifndef PYTHONIC_NUMPY_RANDOM_F_HPP
#define PYTHONIC_NUMPY_RANDOM_F_HPP

#include "pythonic/include/numpy/random/f.hpp"
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
    types::ndarray<double, pS> f(double dfnum, double dfden, pS const &shape)
    {
      types::ndarray<double, pS> result{shape, types::none_type()};
      std::chi_squared_distribution<double> distribution{dfnum};
      std::chi_squared_distribution<double> distribution2{dfden};
      std::generate(result.fbegin(), result.fend(), [&]() {
        return (distribution(details::generator) * dfden) /
               (distribution2(details::generator) * dfnum);
      });
      return result;
    }

    inline auto f(double dfnum, double dfden, long size)
        -> decltype(f(dfnum, dfden, types::array_tuple<long, 1>{{size}}))
    {
      return f(dfnum, dfden, types::array_tuple<long, 1>{{size}});
    }

    inline double f(double dfnum, double dfden, types::none_type d)
    {
      return (std::chi_squared_distribution<double>{dfnum}(details::generator)*dfden) /
             (std::chi_squared_distribution<double>{dfden}(details::generator)*dfnum);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
