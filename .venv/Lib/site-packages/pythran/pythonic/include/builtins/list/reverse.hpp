#ifndef PYTHONIC_INCLUDE_BUILTIN_LIST_REVERSE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_LIST_REVERSE_HPP

#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace list
  {

    template <class T>
    types::none_type reverse(types::list<T> &seq);

    DEFINE_FUNCTOR(pythonic::builtins::list, reverse);
  } // namespace list
} // namespace builtins
PYTHONIC_NS_END
#endif
