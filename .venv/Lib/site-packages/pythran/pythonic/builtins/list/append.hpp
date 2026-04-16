#ifndef PYTHONIC_BUILTIN_LIST_APPEND_HPP
#define PYTHONIC_BUILTIN_LIST_APPEND_HPP

#include "pythonic/include/builtins/list/append.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/types/NoneType.hpp"
#include "pythonic/types/list.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace list
  {

    template <class T, class F>
    types::none_type append(types::list<T> &seq, F &&value)
    {
      seq.push_back(std::forward<F>(value));
      return builtins::None;
    }
    template <class T, class F>
    types::none_type append(types::list<T> &&seq, F &&value)
    {
      seq.push_back(std::forward<F>(value));
      return builtins::None;
    }
  } // namespace list
} // namespace builtins
PYTHONIC_NS_END
#endif
