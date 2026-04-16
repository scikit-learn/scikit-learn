#ifndef PYTHONIC_INCLUDE_BUILTIN_SET_SYMMETRICDIFFERENCEUPDATE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SET_SYMMETRICDIFFERENCEUPDATE_HPP

#include "pythonic/include/types/set.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {

    template <typename T, typename U>
    types::none_type symmetric_difference_update(types::set<T> &set, U const &other);

    template <typename T, typename U>
    types::none_type symmetric_difference_update(types::set<T> &&set, U const &other);

    template <typename U>
    types::none_type symmetric_difference_update(types::empty_set const &set, U const &other);

    DEFINE_FUNCTOR(pythonic::builtins::set, symmetric_difference_update);
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
