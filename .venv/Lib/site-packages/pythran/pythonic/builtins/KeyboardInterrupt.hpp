#ifndef PYTHONIC_BUILTIN_KEYBOARDINTERRUPT_HPP
#define PYTHONIC_BUILTIN_KEYBOARDINTERRUPT_HPP

#include "pythonic/include/builtins/KeyboardInterrupt.hpp"

#include "pythonic/types/exceptions.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  PYTHONIC_EXCEPTION_IMPL(KeyboardInterrupt)
}
PYTHONIC_NS_END

#endif
