#ifndef PYTHONIC_BUILTIN_SET_DIFFERENCEUPDATE_HPP
#define PYTHONIC_BUILTIN_SET_DIFFERENCEUPDATE_HPP

#include "pythonic/include/builtins/set/difference_update.hpp"

#include "pythonic/types/set.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {

    template <typename T, typename... Types>
    types::none_type difference_update(types::set<T> &set, Types const &...others)
    {
      set.difference_update(others...);
      return {};
    }

    template <typename T, typename... Types>
    types::none_type difference_update(types::set<T> &&set, Types const &...others)
    {
      // nothing to be done as we work on rvalue
      return {};
    }

    template <typename... Types>
    types::none_type difference_update(types::empty_set const &set, Types const &...others)
    {
      // nothing can be removed in set
      return {};
    }
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
