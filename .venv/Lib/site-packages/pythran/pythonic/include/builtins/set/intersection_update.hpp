#ifndef PYTHONIC_INCLUDE_BUILTIN_SET_INTERSECTIONUPDATE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SET_INTERSECTIONUPDATE_HPP

#include "pythonic/include/types/set.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {

    template <typename T, typename... Types>
    types::none_type intersection_update(types::set<T> &set, Types const &...others);

    template <typename T, typename... Types>
    types::none_type intersection_update(types::set<T> &&set, Types const &...others);

    template <typename... Types>
    types::none_type intersection_update(types::empty_set &&set, Types const &...others);

    DEFINE_FUNCTOR(pythonic::builtins::set, intersection_update);
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
