#ifndef PYTHONIC_BUILTIN_STR_CAPITALIZE_HPP
#define PYTHONIC_BUILTIN_STR_CAPITALIZE_HPP

#include "pythonic/include/builtins/str/capitalize.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::str capitalize(types::str const &s)
    {
      if (s.empty())
        return s;
      else {
        types::str copy = s;
        copy.chars()[0] = ::toupper(s.chars()[0]);
        std::transform(s.chars().begin() + 1, s.chars().end(), copy.chars().begin() + 1, ::tolower);
        return copy;
      }
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
