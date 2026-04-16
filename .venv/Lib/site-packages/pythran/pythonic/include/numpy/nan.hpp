#ifndef PYTHONIC_INCLUDE_NUMPY_NAN_HPP
#define PYTHONIC_INCLUDE_NUMPY_NAN_HPP

#include <limits>

PYTHONIC_NS_BEGIN

namespace numpy
{
  double const nan = std::numeric_limits<double>::quiet_NaN();
}
PYTHONIC_NS_END

#endif
