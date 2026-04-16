#ifndef PYTHONIC_INCLUDE_RANDOM_SEED_HPP
#define PYTHONIC_INCLUDE_RANDOM_SEED_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/random/random.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace random
{
  types::none_type seed(long s);
  types::none_type seed();

  DEFINE_FUNCTOR(pythonic::random, seed);
} // namespace random

PYTHONIC_NS_END

#endif
