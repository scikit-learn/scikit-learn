#ifndef PYTHONIC_BUILTIN_SET_SYMMETRICDIFFERENCEUPDATE_HPP
#define PYTHONIC_BUILTIN_SET_SYMMETRICDIFFERENCEUPDATE_HPP

#include "pythonic/include/builtins/set/symmetric_difference_update.hpp"

#include "pythonic/types/set.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {

    template <typename T, typename U>
    types::none_type symmetric_difference_update(types::set<T> &set, U const &other)
    {
      set.symmetric_difference_update(other);
      return {};
    }

    template <typename T, typename U>
    types::none_type symmetric_difference_update(types::set<T> &&set, U const &other)
    {
      // nothing to be done on rvalue
      return {};
    }

    template <typename U>
    types::none_type symmetric_difference_update(types::empty_set const &set, U const &other)
    {
      // nothing otherwise empty_set have ! its correct type.
      return {};
    }
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
