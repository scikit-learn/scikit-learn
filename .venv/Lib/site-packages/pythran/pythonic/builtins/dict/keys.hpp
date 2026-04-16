#ifndef PYTHONIC_BUILTIN_DICT_KEYS_HPP
#define PYTHONIC_BUILTIN_DICT_KEYS_HPP

#include "pythonic/include/builtins/dict/keys.hpp"

#include "pythonic/types/dict.hpp"
#include "pythonic/types/list.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace dict
  {

    // We need a copy here for lvalue like :
    // for i in {"a": "b", "c": "d"}.keys():
    //     pass
    template <class D>
    auto keys(D &&d) -> decltype(std::forward<D>(d).keys())
    {
      return std::forward<D>(d).keys();
    }
  } // namespace dict
} // namespace builtins
PYTHONIC_NS_END

#endif
