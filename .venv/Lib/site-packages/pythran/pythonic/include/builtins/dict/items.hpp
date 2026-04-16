#ifndef PYTHONIC_INCLUDE_BUILTIN_DICT_ITEMS_HPP
#define PYTHONIC_INCLUDE_BUILTIN_DICT_ITEMS_HPP

#include "pythonic/include/types/dict.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <tuple>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace dict
  {

    template <class D>
    auto items(D &&d) -> decltype(std::forward<D>(d).items());

    DEFINE_FUNCTOR(pythonic::builtins::dict, items);
  } // namespace dict
} // namespace builtins
PYTHONIC_NS_END

#endif
