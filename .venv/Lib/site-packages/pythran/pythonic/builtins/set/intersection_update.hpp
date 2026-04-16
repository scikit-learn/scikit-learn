#ifndef PYTHONIC_BUILTIN_SET_INTERSECTIONUPDATE_HPP
#define PYTHONIC_BUILTIN_SET_INTERSECTIONUPDATE_HPP

#include "pythonic/include/builtins/set/intersection_update.hpp"

#include "pythonic/types/set.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {

    template <typename T, typename... Types>
    types::none_type intersection_update(types::set<T> &set, Types const &...others)
    {
      set.intersection_update(others...);
      return {};
    }

    template <typename T, typename... Types>
    types::none_type intersection_update(types::set<T> &&set, Types const &...others)
    {
      // If it is an rvalue, we don't really want to update
      return {};
    }

    template <typename... Types>
    types::none_type intersection_update(types::empty_set &&set, Types const &...others)
    {
      // If it is an empty_set, it is ! really updated otherwise we have a
      // typing issue
      return {};
    }
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
