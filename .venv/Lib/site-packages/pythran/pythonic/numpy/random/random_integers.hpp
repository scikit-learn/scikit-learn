#ifndef PYTHONIC_NUMPY_RANDOM_RANDOM_INTEGERS_HPP
#define PYTHONIC_NUMPY_RANDOM_RANDOM_INTEGERS_HPP

#include "pythonic/include/numpy/random/random_integers.hpp"

#include "pythonic/numpy/random/randint.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    template <class T>
    auto random_integers(long min, long max, T &&size)
        -> decltype(randint(min, max, std::forward<T>(size)))
    {
      return randint(min, max + 1, std::forward<T>(size));
    }

    inline long random_integers(long max)
    {
      return randint(1, max + 1);
    }

    inline long random_integers(long min, long max)
    {
      return randint(min, max + 1);
    }
  } // namespace random
} // namespace numpy
PYTHONIC_NS_END

#endif
