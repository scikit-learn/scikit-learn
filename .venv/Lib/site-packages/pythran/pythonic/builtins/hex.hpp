#ifndef PYTHONIC_BUILTIN_HEX_HPP
#define PYTHONIC_BUILTIN_HEX_HPP

#include "pythonic/include/builtins/hex.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

#include <sstream>

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class T>
  types::str hex(T const &v)
  {
    std::ostringstream oss;
    oss << "0x" << std::hex << v;
    return oss.str();
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
