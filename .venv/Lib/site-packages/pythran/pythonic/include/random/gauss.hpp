#ifndef PYTHONIC_INCLUDE_RANDOM_GAUSS_HPP
#define PYTHONIC_INCLUDE_RANDOM_GAUSS_HPP

#include "pythonic/include/random/random.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace random
{

  double gauss(double mu, double sigma);

  DEFINE_FUNCTOR(pythonic::random, gauss);
} // namespace random
PYTHONIC_NS_END

#endif
