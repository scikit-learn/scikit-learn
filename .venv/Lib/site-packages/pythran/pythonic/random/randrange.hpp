#ifndef PYTHONIC_RANDOM_RANDRANGE_HPP
#define PYTHONIC_RANDOM_RANDRANGE_HPP

#include "pythonic/include/random/randrange.hpp"

#include "pythonic/random/random.hpp"
#include "pythonic/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace random
{
  inline long randrange(long stop)
  {
    return long(random() * stop);
  }

  inline long randrange(long start, long stop)
  {
    return start + long(random() * (stop - start));
  }

  inline long randrange(long start, long stop, long step)
  {
    return start + step * long((random() * (stop - start)) / std::abs(step));
  }
} // namespace random
PYTHONIC_NS_END

#endif
