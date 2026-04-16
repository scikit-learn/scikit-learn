#ifndef PYTHONIC_RANDOM_UNIFORM_HPP
#define PYTHONIC_RANDOM_UNIFORM_HPP

#include "pythonic/include/random/uniform.hpp"

#include "pythonic/random/random.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace random
{
  inline double uniform(double a, double b)
  {
    return a + (b - a) * random();
  }
} // namespace random
PYTHONIC_NS_END

#endif
