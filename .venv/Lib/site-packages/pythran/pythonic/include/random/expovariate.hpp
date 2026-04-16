#ifndef PYTHONIC_INCLUDE_RANDOM_EXPOVARIATE_HPP
#define PYTHONIC_INCLUDE_RANDOM_EXPOVARIATE_HPP

#include "pythonic/include/random/random.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace random
{
  double expovariate(double l);

  DEFINE_FUNCTOR(pythonic::random, expovariate);
} // namespace random
PYTHONIC_NS_END

#endif
