#ifndef PYTHONIC_BUILTIN_SET_ADD_HPP
#define PYTHONIC_BUILTIN_SET_ADD_HPP

#include "pythonic/include/builtins/set/add.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/types/NoneType.hpp"
#include "pythonic/types/set.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {

    template <class T, class F>
    types::none_type add(types::set<T> &s, F const &value)
    {
      s.add(value);
      return builtins::None;
    }

    template <class T, class F>
    types::none_type add(types::set<T> &&s, F const &value)
    {
      s.add(value);
      return builtins::None;
    }
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
