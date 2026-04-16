#ifndef PYTHONIC_RANDOM_RANDOM_HPP
#define PYTHONIC_RANDOM_RANDOM_HPP

#include "pythonic/include/random/random.hpp"

#include "pythonic/utils/functor.hpp"
#include <random>

PYTHONIC_NS_BEGIN

namespace random
{
  inline double random()
  {
    static std::uniform_real_distribution<> uniform_distrib(0.0, 1.0);
    return uniform_distrib(__random_generator);
  }
} // namespace random
PYTHONIC_NS_END

#endif
