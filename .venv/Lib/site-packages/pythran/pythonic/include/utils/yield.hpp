#ifndef PYTHRAN_INCLUDE_UTILS_YIELD_HPP
#define PYTHRAN_INCLUDE_UTILS_YIELD_HPP

/*
 * This contains base class for yielders
 */

#include "pythonic/include/types/generator.hpp"

PYTHONIC_NS_BEGIN
class yielder
{
public:
  yielder();

  bool operator!=(yielder const &other) const;
  bool operator==(yielder const &other) const;

  long __generator_state;
};
PYTHONIC_NS_END

#endif
