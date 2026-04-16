#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_FIND_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_FIND_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    long find(types::str const &s, types::str const &value, long start, long end);

    long find(types::str const &s, types::str const &value, long start);

    long find(types::str const &s, types::str const &value);

    DEFINE_FUNCTOR(pythonic::builtins::str, find);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
