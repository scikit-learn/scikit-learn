#ifndef PYTHONIC_INCLUDE_BUILTIN_DICT_KEYS_HPP
#define PYTHONIC_INCLUDE_BUILTIN_DICT_KEYS_HPP

#include "pythonic/include/types/dict.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace dict
  {

    template <class D>
    auto keys(D &&d) -> decltype(std::forward<D>(d).keys());

    DEFINE_FUNCTOR(pythonic::builtins::dict, keys);
  } // namespace dict
} // namespace builtins
PYTHONIC_NS_END

#endif
