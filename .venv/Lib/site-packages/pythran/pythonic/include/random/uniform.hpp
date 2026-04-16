#ifndef PYTHONIC_INCLUDE_RANDOM_UNIFORM_HPP
#define PYTHONIC_INCLUDE_RANDOM_UNIFORM_HPP

#include "pythonic/include/random/random.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace random
{
  double uniform(double a, double b);

  DEFINE_FUNCTOR(pythonic::random, uniform);
} // namespace random
PYTHONIC_NS_END

#endif
