#ifndef PYTHONIC_NUMPY_RANDOM_RANDN_HPP
#define PYTHONIC_NUMPY_RANDOM_RANDN_HPP

#include "pythonic/include/numpy/random/randn.hpp"

#include "pythonic/numpy/random/standard_normal.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {

    template <class... T>
    types::ndarray<double, types::array_tuple<long, sizeof...(T)>> randn(T... shape)
    {
      return standard_normal(types::array_tuple<long, sizeof...(T)>{{shape...}});
    }

    inline double randn()
    {
      return standard_normal();
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
