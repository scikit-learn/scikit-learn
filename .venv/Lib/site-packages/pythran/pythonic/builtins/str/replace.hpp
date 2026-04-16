#ifndef PYTHONIC_BUILTIN_STR_REPLACE_HPP
#define PYTHONIC_BUILTIN_STR_REPLACE_HPP

#include "pythonic/include/builtins/str/replace.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/allocate.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::str replace(types::str const &self, types::str const &old_pattern,
                       types::str const &new_pattern, long count)
    {
      char const *needle = old_pattern.c_str();
      char const *new_needle = new_pattern.c_str();
      char const *new_needle_end = new_needle + new_pattern.size();
      char const *haystack = self.c_str();

      char const *haystack_next = strstr(haystack, needle);
      if (!count || !haystack_next) {
        return {haystack};
      } else {
        size_t n = 1 + std::max(self.size(),
                                self.size() * (1 + new_pattern.size()) / (1 + old_pattern.size()));

        char *buffer = utils::allocate<char>(n);
        char *iter = buffer;
        do {
          iter = std::copy(haystack, haystack_next, iter);
          iter = std::copy(new_needle, new_needle_end, iter);
          --count;
          haystack = haystack_next + old_pattern.size();
          assert(size_t(iter - buffer) < n);
        } while (count && (haystack_next = strstr(haystack, needle)));

        iter = std::copy(haystack, self.c_str() + self.size(), iter);
        types::str replaced(buffer, iter - buffer);
        utils::deallocate(buffer);
        return replaced;
      }
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
