#ifndef PYTHONIC_INCLUDE_BUILTIN_SET_DIFFERENCEUPDATE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SET_DIFFERENCEUPDATE_HPP

#include "pythonic/include/types/set.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {

    template <typename T, typename... Types>
    types::none_type difference_update(types::set<T> &set, Types const &...others);

    template <typename T, typename... Types>
    types::none_type difference_update(types::set<T> &&set, Types const &...others);

    template <typename... Types>
    types::none_type difference_update(types::empty_set const &set, Types const &...others);

    DEFINE_FUNCTOR(pythonic::builtins::set, difference_update);
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
