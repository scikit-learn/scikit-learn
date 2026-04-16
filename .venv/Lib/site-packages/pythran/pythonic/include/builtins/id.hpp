#ifndef PYTHONIC_INCLUDE_BUILTIN_ID_HPP
#define PYTHONIC_INCLUDE_BUILTIN_ID_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class T>
  long id(T const &t);

  long id(long const &t);
  long id(double const &t);
  long id(bool const &t);

  DEFINE_FUNCTOR(pythonic::builtins, id);
} // namespace builtins
PYTHONIC_NS_END

#endif
