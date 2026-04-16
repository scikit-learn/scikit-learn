#ifndef PYTHONIC_INCLUDE_BUILTIN_SET_ISDISJOINT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SET_ISDISJOINT_HPP

#include "pythonic/include/types/set.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {
    template <class T, class U>
    bool isdisjoint(types::set<T> const &calling_set, U const &arg_set);

    template <class U>
    bool isdisjoint(types::empty_set const &calling_set, U const &arg_set);

    DEFINE_FUNCTOR(pythonic::builtins::set, isdisjoint);
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
