#ifndef PYTHONIC_RANDOM_EXPOVARIATE_HPP
#define PYTHONIC_RANDOM_EXPOVARIATE_HPP

#include "pythonic/include/random/expovariate.hpp"

#include "pythonic/random/random.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace random
{
  inline double expovariate(double l)
  {
    return std::exponential_distribution<>(l)(__random_generator);
  }
} // namespace random
PYTHONIC_NS_END

#endif
