#ifndef PYTHONIC_INCLUDE_DISPATCH_INSERT_HPP
#define PYTHONIC_INCLUDE_DISPATCH_INSERT_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{

  template <class Any, class Arg>
  types::none_type insert(Any &&any, long index, Arg &&arg0);

  DEFINE_FUNCTOR(pythonic::__dispatch__, insert);
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
