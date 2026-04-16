#ifndef PYTHONIC_NUMPY_RANDOM_STANDARD_GAMMA_HPP
#define PYTHONIC_NUMPY_RANDOM_STANDARD_GAMMA_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/standard_gamma.hpp"

#include "pythonic/numpy/random/gamma.hpp"
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
    types::ndarray<double, pS> standard_gamma(double s, pS const &shape)
    {
      return gamma(s, 1., shape);
    }

    inline auto standard_gamma(double s, long size)
        -> decltype(standard_gamma(s, types::array_tuple<long, 1>{{size}}))
    {
      return standard_gamma(s, types::array_tuple<long, 1>{{size}});
    }

    inline double standard_gamma(double s, types::none_type d)
    {
      return gamma(s, 1., d);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
