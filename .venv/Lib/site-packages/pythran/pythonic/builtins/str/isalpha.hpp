#ifndef PYTHONIC_BUILTIN_STR_ISALPHA_HPP
#define PYTHONIC_BUILTIN_STR_ISALPHA_HPP

#include "pythonic/include/builtins/str/isalpha.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    bool isalpha(types::str const &s)
    {
      return !s.empty() &&
             std::all_of(s.chars().begin(), s.chars().end(), (int (*)(int))std::isalpha);
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
