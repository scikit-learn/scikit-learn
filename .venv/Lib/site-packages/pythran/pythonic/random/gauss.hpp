#ifndef PYTHONIC_RANDOM_GAUSS_HPP
#define PYTHONIC_RANDOM_GAUSS_HPP

#include "pythonic/include/random/gauss.hpp"

#include "pythonic/random/random.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace random
{

  inline double gauss(double mu, double sigma)
  {
    return std::normal_distribution<>(mu, sigma)(__random_generator);
  }
} // namespace random
PYTHONIC_NS_END

#endif
