#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_RANDN_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_RANDN_HPP

#include "pythonic/include/numpy/random/generator.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    template <class... T>
    types::ndarray<double, types::array_tuple<long, sizeof...(T)>> randn(T... shape);

    double randn();

    DEFINE_FUNCTOR(pythonic::numpy::random, randn);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
