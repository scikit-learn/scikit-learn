#ifndef PYTHONIC_INCLUDE_TIME_SLEEP_HPP
#define PYTHONIC_INCLUDE_TIME_SLEEP_HPP

#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace time
{
  types::none_type sleep(double const value);

  DEFINE_FUNCTOR(pythonic::time, sleep)
} // namespace time
PYTHONIC_NS_END

#endif
