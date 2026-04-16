#ifndef PYTHONIC_BUILTIN_STR_FIND_HPP
#define PYTHONIC_BUILTIN_STR_FIND_HPP

#include "pythonic/include/builtins/str/find.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    inline long find(types::str const &s, types::str const &value, long start, long end)
    {
      if (end < 0)
        end += s.size();
      long a = s.find(value, start);
      return (a > end) ? -1 : a;
    }

    inline long find(types::str const &s, types::str const &value, long start)
    {
      return find(s, value, start, s.size());
    }

    inline long find(types::str const &s, types::str const &value)
    {
      return find(s, value, 0, s.size());
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
