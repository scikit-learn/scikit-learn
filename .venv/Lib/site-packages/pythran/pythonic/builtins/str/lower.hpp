#ifndef PYTHONIC_BUILTIN_STR_LOWER_HPP
#define PYTHONIC_BUILTIN_STR_LOWER_HPP

#include "pythonic/include/builtins/str/lower.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::str lower(types::str const &s)
    {
      types::str copy = s;
      std::transform(s.chars().begin(), s.chars().end(), copy.chars().begin(), ::tolower);
      return copy;
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
