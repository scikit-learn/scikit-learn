#ifndef PYTHONIC_BUILTIN_STR_RSTRIP_HPP
#define PYTHONIC_BUILTIN_STR_RSTRIP_HPP

#include "pythonic/include/builtins/str/rstrip.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::str rstrip(types::str const &self, types::str const &to_del)
    {
      auto chars = self.chars();
      auto stop = self.find_last_not_of(to_del);
      if (stop < 0)
        return {};
      return {chars.begin(), chars.begin() + stop + 1};
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
