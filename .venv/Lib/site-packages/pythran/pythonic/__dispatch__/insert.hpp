#ifndef PYTHONIC_DISPATCH_INSERT_HPP
#define PYTHONIC_DISPATCH_INSERT_HPP

#include "pythonic/include/__dispatch__/insert.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{

  template <class Any, class Arg>
  types::none_type insert(Any &&any, long index, Arg &&arg)
  {
    index = index % (1 + any.size()); // +1 because we want to be able to insert
                                      // at the end of any
    if (index < 0)
      index += any.size();
    any.insert(index, std::forward<Arg>(arg));
    return builtins::None;
  }
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
