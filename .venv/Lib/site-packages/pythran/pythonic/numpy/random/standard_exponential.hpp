#ifndef PYTHONIC_NUMPY_RANDOM_STANDARD_EXPONENTIAL_HPP
#define PYTHONIC_NUMPY_RANDOM_STANDARD_EXPONENTIAL_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/standard_exponential.hpp"

#include "pythonic/numpy/random/exponential.hpp"
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
    types::ndarray<double, pS> standard_exponential(pS const &shape)
    {
      return exponential(1., shape);
    }

    inline auto standard_exponential(long size)
        -> decltype(standard_exponential(types::array_tuple<long, 1>{{size}}))
    {
      return standard_exponential(types::array_tuple<long, 1>{{size}});
    }

    inline double standard_exponential(types::none_type d)
    {
      return exponential(1., d);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
