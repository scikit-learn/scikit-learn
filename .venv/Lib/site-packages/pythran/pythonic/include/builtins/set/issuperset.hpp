#ifndef PYTHONIC_INCLUDE_BUILTIN_SET_ISSUPERSET_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SET_ISSUPERSET_HPP

#include "pythonic/include/types/set.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {

    template <class T, class U>
    bool issuperset(types::set<T> const &set, U const &other);

    template <class U>
    bool issuperset(types::empty_set const &set, U const &other);

    DEFINE_FUNCTOR(pythonic::builtins::set, issuperset);
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
