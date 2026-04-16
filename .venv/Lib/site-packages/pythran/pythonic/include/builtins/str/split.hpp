#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_SPLIT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_SPLIT_HPP

#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::list<types::str> split(types::str const &in, types::str const &sep, long maxsplit = -1);

    types::list<types::str> split(types::str const &s, types::none_type const & = {},
                                  long maxsplit = -1);

    DEFINE_FUNCTOR(pythonic::builtins::str, split);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
