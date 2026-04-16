#ifndef PYTHONIC_BUILTIN_STR_LSTRIP_HPP
#define PYTHONIC_BUILTIN_STR_LSTRIP_HPP

#include "pythonic/include/builtins/str/lstrip.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::str lstrip(types::str const &self, types::str const &to_del)
    {
      auto chars = self.chars();
      auto stop = self.find_first_not_of(to_del);
      if (stop < 0)
        return {};
      else
        return {chars.begin() + stop, chars.end()};
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
