#ifndef PYTHONIC_INCLUDE_BUILTIN_LIST_APPEND_HPP
#define PYTHONIC_INCLUDE_BUILTIN_LIST_APPEND_HPP

#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace list
  {

    template <class T, class F>
    types::none_type append(types::list<T> &seq, F &&value);

    template <class T, class F>
    types::none_type append(types::list<T> &&seq, F &&value);

    template <class F>
    types::none_type append(types::empty_list &seq, F &&value);

    DEFINE_FUNCTOR(pythonic::builtins::list, append);
  } // namespace list
} // namespace builtins
PYTHONIC_NS_END
#endif
