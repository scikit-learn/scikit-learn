#ifndef PYTHONIC_INCLUDE_BUILTIN_DIVMOD_HPP
#define PYTHONIC_INCLUDE_BUILTIN_DIVMOD_HPP

#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class T0, class T1>
  auto divmod(T0 const &t0, T1 const &t1) // other types are left over
      -> decltype(types::make_tuple(t0 / t1, t0 % t1));

  DEFINE_FUNCTOR(pythonic::builtins, divmod);
} // namespace builtins
PYTHONIC_NS_END

#endif
