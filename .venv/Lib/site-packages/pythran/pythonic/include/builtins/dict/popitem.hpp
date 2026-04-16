#ifndef PYTHONIC_INCLUDE_BUILTIN_DICT_POPITEM_HPP
#define PYTHONIC_INCLUDE_BUILTIN_DICT_POPITEM_HPP

#include "pythonic/include/types/dict.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <tuple>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace dict
  {

    template <class D>
    auto popitem(D &&d) -> decltype(std::forward<D>(d).popitem());

    DEFINE_FUNCTOR(pythonic::builtins::dict, popitem);
  } // namespace dict
} // namespace builtins
PYTHONIC_NS_END

#endif
