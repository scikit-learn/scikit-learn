#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_RAND_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_RAND_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    template <class... T>
    types::ndarray<double, types::array_tuple<long, sizeof...(T)>> rand(T... shape);

    double rand();

    DEFINE_FUNCTOR(pythonic::numpy::random, rand);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
