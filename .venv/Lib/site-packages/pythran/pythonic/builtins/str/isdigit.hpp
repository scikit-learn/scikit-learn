#ifndef PYTHONIC_BUILTIN_STR_ISDIGIT_HPP
#define PYTHONIC_BUILTIN_STR_ISDIGIT_HPP

#include "pythonic/include/builtins/str/isdigit.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

#include <cctype>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    bool isdigit(types::str const &s)
    {
      return !s.empty() &&
             std::all_of(s.chars().begin(), s.chars().end(), (int (*)(int))std::isdigit);
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
