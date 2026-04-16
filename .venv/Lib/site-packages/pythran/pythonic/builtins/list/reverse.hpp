#ifndef PYTHONIC_BUILTIN_LIST_REVERSE_HPP
#define PYTHONIC_BUILTIN_LIST_REVERSE_HPP

#include "pythonic/include/builtins/list/reverse.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/types/NoneType.hpp"
#include "pythonic/types/list.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace list
  {

    template <class T>
    types::none_type reverse(types::list<T> &seq)
    {
      std::reverse(seq.begin(), seq.end());
      return builtins::None;
    }
  } // namespace list
} // namespace builtins
PYTHONIC_NS_END
#endif
