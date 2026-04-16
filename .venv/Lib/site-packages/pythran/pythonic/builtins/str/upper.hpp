#ifndef PYTHONIC_BUILTIN_STR_UPPER_HPP
#define PYTHONIC_BUILTIN_STR_UPPER_HPP

#include "pythonic/include/builtins/str/upper.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {

    types::str upper(types::str const &s)
    {
      types::str copy = s;
      std::transform(s.chars().begin(), s.chars().end(), copy.chars().begin(), ::toupper);
      return copy;
    }
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
