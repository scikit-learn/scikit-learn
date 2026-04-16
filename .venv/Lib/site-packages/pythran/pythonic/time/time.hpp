#ifndef PYTHONIC_TIME_TIME_HPP
#define PYTHONIC_TIME_TIME_HPP

#include "pythonic/include/time/time.hpp"
#include "pythonic/utils/functor.hpp"

#include <chrono>

PYTHONIC_NS_BEGIN

namespace time
{

  double time()
  {
    std::chrono::time_point<std::chrono::steady_clock> tp = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch()).count() /
           1000.;
  }
} // namespace time
PYTHONIC_NS_END

#endif
