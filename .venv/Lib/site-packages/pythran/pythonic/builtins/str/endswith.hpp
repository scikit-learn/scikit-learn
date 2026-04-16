#ifndef PYTHONIC_BUILTIN_STR_ENDSWITH_HPP
#define PYTHONIC_BUILTIN_STR_ENDSWITH_HPP

#include "pythonic/include/builtins/str/endswith.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    bool endswith(types::str const &s, types::str const &suffix, long start, long end)
    {
      if (end == -1)
        end = s.size();
      long rstart = end - suffix.size();
      return rstart >= start && s.compare(rstart, suffix.size(), suffix) == 0;
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
