#ifndef PYTHONIC_INCLUDE_MATH_FMOD_HPP
#define PYTHONIC_INCLUDE_MATH_FMOD_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(fmod, std::fmod);
}
PYTHONIC_NS_END

#endif
