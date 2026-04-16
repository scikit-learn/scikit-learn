#ifndef PYTHONIC_NUMPY_RANDOM_STANDARD_NORMAL_HPP
#define PYTHONIC_NUMPY_RANDOM_STANDARD_NORMAL_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/numpy/random/standard_normal.hpp"

#include "pythonic/numpy/random/normal.hpp"
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
    types::ndarray<double, pS> standard_normal(pS const &shape)
    {
      return normal(0., 1., shape);
    }

    inline auto standard_normal(long size)
        -> decltype(standard_normal(types::array_tuple<long, 1>{{size}}))
    {
      return standard_normal(types::array_tuple<long, 1>{{size}});
    }

    inline double standard_normal(types::none_type d)
    {
      return normal(0., 1., d);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
