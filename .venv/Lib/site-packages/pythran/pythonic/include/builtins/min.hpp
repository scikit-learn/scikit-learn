#ifndef PYTHONIC_INCLUDE_BUILTIN_MIN_HPP
#define PYTHONIC_INCLUDE_BUILTIN_MIN_HPP

#include "pythonic/include/builtins/minmax.hpp"
#include "pythonic/include/operator_/gt.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class... Types>
  auto min(Types &&...values)
      -> decltype(details::minmax(operator_::functor::gt{}, std::forward<Types>(values)...));

  DEFINE_FUNCTOR(pythonic::builtins, min);
} // namespace builtins
PYTHONIC_NS_END

#endif
