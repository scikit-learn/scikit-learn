#ifndef PYTHONIC_BUILTIN_STR_STRIP_HPP
#define PYTHONIC_BUILTIN_STR_STRIP_HPP

#include "pythonic/include/builtins/str/strip.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {
    types::str strip(types::str const &self, types::str const &to_del)
    {
      if (!self)
        return self;
      auto first = self.find_first_not_of(to_del);
      if (first == -1)
        return types::str();
      else
        return types::str(self.chars().begin() + first,
                          self.chars().begin() + self.find_last_not_of(to_del) + 1);
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
