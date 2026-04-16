#ifndef PYTHONIC_INCLUDE_MATH_HYPOT_HPP
#define PYTHONIC_INCLUDE_MATH_HYPOT_HPP

#include "pythonic/include/utils/functor.hpp"
#include <cmath>

#undef hypot
// This is a windows defined macro that clash with std::hypot && our hypot
// function

PYTHONIC_NS_BEGIN

namespace math
{
  DEFINE_FUNCTOR_2(hypot, std::hypot);
}
PYTHONIC_NS_END

#endif
