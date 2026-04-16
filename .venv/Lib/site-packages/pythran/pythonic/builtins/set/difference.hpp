#ifndef PYTHONIC_BUILTIN_SET_DIFFERENCE_HPP
#define PYTHONIC_BUILTIN_SET_DIFFERENCE_HPP

#include "pythonic/include/builtins/set/difference.hpp"

#include "pythonic/types/set.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {

    template <typename T, typename... Types>
    types::set<T> difference(types::set<T> const &set, Types const &...others)
    {
      return set.difference(others...);
    }

    template <typename T, typename... Types>
    types::set<T> difference(types::set<T> &&set, Types const &...others)
    {
      set.difference_update(others...);
      return set;
    }

    template <typename... Types>
    types::empty_set difference(types::empty_set const &set, Types const &...others)
    {
      return types::empty_set();
    }

    template <typename T>
    types::set<T> difference(types::set<T> const &set)
    {
      return set;
    }

    template <typename T>
    types::set<T> difference(types::set<T> &&set)
    {
      return set;
    }

    inline types::empty_set difference(types::empty_set const &set)
    {
      return types::empty_set();
    }
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
