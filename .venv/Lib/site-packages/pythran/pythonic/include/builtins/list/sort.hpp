#ifndef PYTHONIC_INCLUDE_BUILTIN_LIST_SORT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_LIST_SORT_HPP

#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace list
  {

    template <class T>
    types::none_type sort(types::list<T> &seq);

    template <class T, class K>
    types::none_type sort(types::list<T> &seq, K key);

    DEFINE_FUNCTOR(pythonic::builtins::list, sort);
  } // namespace list
} // namespace builtins
PYTHONIC_NS_END
#endif
