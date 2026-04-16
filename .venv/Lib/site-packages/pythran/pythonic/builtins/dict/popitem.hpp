#ifndef PYTHONIC_BUILTIN_DICT_POPITEM_HPP
#define PYTHONIC_BUILTIN_DICT_POPITEM_HPP

#include "pythonic/include/builtins/dict/popitem.hpp"

#include "pythonic/types/dict.hpp"
#include "pythonic/utils/functor.hpp"

#include <tuple>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace dict
  {

    template <class D>
    auto popitem(D &&d) -> decltype(std::forward<D>(d).popitem())
    {
      return std::forward<D>(d).popitem();
    }
  } // namespace dict
} // namespace builtins
PYTHONIC_NS_END

#endif
