#ifndef PYTHONIC_INCLUDE_NUMPY_NINF_HPP
#define PYTHONIC_INCLUDE_NUMPY_NINF_HPP

#include <limits>

PYTHONIC_NS_BEGIN

namespace numpy
{
  double const NINF = -std::numeric_limits<double>::infinity();
}
PYTHONIC_NS_END

#endif
