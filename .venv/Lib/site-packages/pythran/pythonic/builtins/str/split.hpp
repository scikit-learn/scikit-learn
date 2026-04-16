#ifndef PYTHONIC_BUILTIN_STR_SPLIT_HPP
#define PYTHONIC_BUILTIN_STR_SPLIT_HPP

#include "pythonic/include/builtins/str/split.hpp"

#include "pythonic/builtins/str/strip.hpp"
#include "pythonic/types/NoneType.hpp"
#include "pythonic/types/list.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    inline types::list<types::str> split(types::str const &in, types::str const &sep, long maxsplit)
    {
      types::str s = strip(in);
      types::list<types::str> res(0);
      if (s.empty())
        return res;

      size_t current = 0;
      size_t next = 0;
      long numsplit = 0;
      while (next != types::str::npos && (numsplit++ < maxsplit || maxsplit == -1)) {
        next = s.find_first_of(sep, current);
        res.push_back(s.substr(current, next - current));
        current = next + 1;
      }
      if (next != types::str::npos) {
        current = next + 1;
        res.push_back(s.substr(current, s.size() - current));
      }
      return res;
    }

    inline types::list<types::str> split(types::str const &in, types::none_type const &,
                                         long maxsplit)
    {
      types::str s = strip(in);
      types::list<types::str> res(0);
      if (s.empty())
        return res;

      size_t current = 0;
      size_t next = 0;
      long numsplit = 0;
      while (next != types::str::npos && (numsplit++ < maxsplit || maxsplit == -1)) {
        next = s.find_first_of(" \n\r\t", current);
        // from the pydoc, we skip any blank list
        size_t end = s.find_first_not_of(" \n\r\t", next);
        res.push_back(s.substr(current, next - current));
        current = end;
      }
      if (next != types::str::npos) {
        current = next + 1;
        res.push_back(s.substr(current, s.size() - current));
      }
      return res;
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
