#ifndef PYTHONIC_BUILTIN_FLOAT_ISINTEGER_HPP
#define PYTHONIC_BUILTIN_FLOAT_ISINTEGER_HPP

#include "pythonic/include/builtins/float_/is_integer.hpp"

#include "pythonic/utils/functor.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace float_
  {

    bool is_integer(double d)
    {
      return std::trunc(d) == d;
    }
  } // namespace float_
} // namespace builtins
PYTHONIC_NS_END

#endif
