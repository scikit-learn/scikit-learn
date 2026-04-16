#ifndef PYTHONIC_INCLUDE_RANDOM_RANDINT_HPP
#define PYTHONIC_INCLUDE_RANDOM_RANDINT_HPP

#include "pythonic/include/random/randrange.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace random
{

  long randint(long a, long b);

  DEFINE_FUNCTOR(pythonic::random, randint);
} // namespace random
PYTHONIC_NS_END

#endif
