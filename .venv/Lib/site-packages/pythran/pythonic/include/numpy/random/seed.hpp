#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_SEED_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_SEED_HPP

#include "pythonic/include/numpy/random/generator.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace random
  {
    types::none_type seed(long s);
    types::none_type seed(types::none_type _ = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, seed);
  } // namespace random
} // namespace numpy

PYTHONIC_NS_END

#endif
