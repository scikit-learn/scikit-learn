#ifndef PYTHONIC_INCLUDE_BUILTIN_DICT_VALUES_HPP
#define PYTHONIC_INCLUDE_BUILTIN_DICT_VALUES_HPP

#include "pythonic/include/types/dict.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace dict
  {
    template <class D>
    auto values(D &&d) -> decltype(std::forward<D>(d).values());

    DEFINE_FUNCTOR(pythonic::builtins::dict, values);
  } // namespace dict
} // namespace builtins
PYTHONIC_NS_END

#endif
