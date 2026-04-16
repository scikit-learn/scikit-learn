#ifndef PYTHONIC_BUILTIN_DICT_VALUES_HPP
#define PYTHONIC_BUILTIN_DICT_VALUES_HPP

#include "pythonic/include/builtins/dict/values.hpp"

#include "pythonic/types/dict.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace dict
  {

    template <class D>
    auto values(D &&d) -> decltype(std::forward<D>(d).values())
    {
      return std::forward<D>(d).values();
    }
  } // namespace dict
} // namespace builtins
PYTHONIC_NS_END

#endif
