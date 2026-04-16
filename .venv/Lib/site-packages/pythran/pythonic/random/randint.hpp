#ifndef PYTHONIC_RANDOM_RANDINT_HPP
#define PYTHONIC_RANDOM_RANDINT_HPP

#include "pythonic/include/random/randint.hpp"

#include "pythonic/random/randrange.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace random
{

  inline long randint(long a, long b)
  {
    // TODO: It should be implemented with an uniform_int_distribution
    return randrange(a, b + 1);
  }
} // namespace random
PYTHONIC_NS_END

#endif
