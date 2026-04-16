#ifndef PYTHONIC_BUILTIN_MIN_HPP
#define PYTHONIC_BUILTIN_MIN_HPP

#include "pythonic/builtins/minmax.hpp"
#include "pythonic/include/builtins/min.hpp"

#include "pythonic/operator_/gt.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class... Types>
  auto min(Types &&...values)
      -> decltype(details::minmax(operator_::functor::gt{}, std::forward<Types>(values)...))
  {
    return details::minmax(operator_::functor::gt{}, std::forward<Types>(values)...);
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
