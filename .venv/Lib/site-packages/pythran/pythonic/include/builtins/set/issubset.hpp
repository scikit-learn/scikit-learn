#ifndef PYTHONIC_INCLUDE_BUILTIN_SET_ISSUBSET_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SET_ISSUBSET_HPP

#include "pythonic/include/types/set.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {

    template <class T, class U>
    bool issubset(types::set<T> const &set, U const &other);

    template <class U>
    bool issubset(types::empty_set const &set, U const &other);

    DEFINE_FUNCTOR(pythonic::builtins::set, issubset);
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
