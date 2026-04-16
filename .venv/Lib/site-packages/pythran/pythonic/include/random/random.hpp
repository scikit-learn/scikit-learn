#ifndef PYTHONIC_INCLUDE_RANDOM_RANDOM_HPP
#define PYTHONIC_INCLUDE_RANDOM_RANDOM_HPP

#include "pythonic/include/utils/functor.hpp"
#include <random>

PYTHONIC_NS_BEGIN

namespace random
{

  static std::mt19937 __random_generator;

  double random();

  DEFINE_FUNCTOR(pythonic::random, random);
} // namespace random
PYTHONIC_NS_END

#endif
