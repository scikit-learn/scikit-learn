#ifndef PYTHONIC_NUMPY_RANDOM_SEED_HPP
#define PYTHONIC_NUMPY_RANDOM_SEED_HPP

#include "pythonic/builtins/None.hpp"
#include "pythonic/include/numpy/random/seed.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace random
  {

    inline types::none_type seed(long s)
    {
      details::generator.seed(s);
      return builtins::None;
    }

    inline types::none_type seed(types::none_type)
    {
      details::generator.seed();
      return builtins::None;
    }
  } // namespace random
} // namespace numpy

PYTHONIC_NS_END

#endif
