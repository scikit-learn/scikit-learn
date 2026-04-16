#ifndef PYTHRAN_UTILS_YIELD_HPP
#define PYTHRAN_UTILS_YIELD_HPP

#include "pythonic/include/utils/yield.hpp"

/*
 * This contains base class for yielders
 */

#include "pythonic/types/generator.hpp"

PYTHONIC_NS_BEGIN
inline yielder::yielder() : __generator_state(0)
{
}

inline bool yielder::operator!=(yielder const &other) const
{
  return __generator_state != other.__generator_state;
}

inline bool yielder::operator==(yielder const &other) const
{
  return __generator_state == other.__generator_state;
}
PYTHONIC_NS_END

#endif
