#ifndef PYTHONIC_BUILTIN_ORD_HPP
#define PYTHONIC_BUILTIN_ORD_HPP

#include "pythonic/include/builtins/ord.hpp"

#include "pythonic/types/exceptions.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  long ord(types::str const &v)
  {
    if (v.size() != 1)
      throw types::TypeError("ord() expected a character, but string of length " +
                             std::to_string(v.size()) + " found");
    return (long)v.chars()[0];
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
