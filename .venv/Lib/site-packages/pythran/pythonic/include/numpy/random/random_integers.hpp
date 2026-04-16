#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_RANDOM_INTEGERS_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_RANDOM_INTEGERS_HPP

#include "pythonic/include/numpy/random/randint.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <utility>

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    template <class T>
    auto random_integers(long min, long max, T &&size)
        -> decltype(randint(min, max, std::forward<T>(size)));

    long random_integers(long max);

    long random_integers(long min, long max);

    DEFINE_FUNCTOR(pythonic::numpy::random, random_integers);
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
