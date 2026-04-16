#ifndef PYTHONIC_INCLUDE_TIME_TIME_HPP
#define PYTHONIC_INCLUDE_TIME_TIME_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace time
{

  double time();

  DEFINE_FUNCTOR(pythonic::time, time)
} // namespace time
PYTHONIC_NS_END

#endif
