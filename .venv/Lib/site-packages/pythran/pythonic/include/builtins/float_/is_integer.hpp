#ifndef PYTHONIC_INCLUDE_BUILTIN_FLOAT_ISINTEGER_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FLOAT_ISINTEGER_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace float_
  {

    bool is_integer(double d);

    DEFINE_FUNCTOR(pythonic::builtins::float_, is_integer);
  } // namespace float_
} // namespace builtins
PYTHONIC_NS_END

#endif
