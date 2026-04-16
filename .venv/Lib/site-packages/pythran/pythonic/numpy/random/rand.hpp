#ifndef PYTHONIC_NUMPY_RANDOM_RAND_HPP
#define PYTHONIC_NUMPY_RANDOM_RAND_HPP

#include "pythonic/include/numpy/random/rand.hpp"

#include "pythonic/numpy/random/random.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {

    template <class... T>
    types::ndarray<double, types::array_tuple<long, sizeof...(T)>> rand(T... shape)
    {
      return random(types::array_tuple<long, sizeof...(T)>{{shape...}});
    }

    inline double rand()
    {
      return random();
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
