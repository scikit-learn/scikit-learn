#ifndef PYTHONIC_INCLUDE_NUMPY_INF_HPP
#define PYTHONIC_INCLUDE_NUMPY_INF_HPP
#include <limits>

PYTHONIC_NS_BEGIN

namespace numpy
{
  double const inf = std::numeric_limits<double>::infinity();
}
PYTHONIC_NS_END

#endif
