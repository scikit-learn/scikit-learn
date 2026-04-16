#ifndef PYTHONIC_BUILTIN_MAX_HPP
#define PYTHONIC_BUILTIN_MAX_HPP

#include "pythonic/builtins/minmax.hpp"
#include "pythonic/include/builtins/max.hpp"

#include "pythonic/operator_/lt.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class... Types>
  auto max(Types &&...values)
      -> decltype(details::minmax(operator_::functor::lt{}, std::forward<Types>(values)...))
  {
    return details::minmax(operator_::functor::lt{}, std::forward<Types>(values)...);
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
