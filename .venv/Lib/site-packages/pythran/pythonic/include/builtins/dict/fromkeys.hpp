#ifndef PYTHONIC_INCLUDE_BUILTIN_DICT_FROMKEYS_HPP
#define PYTHONIC_INCLUDE_BUILTIN_DICT_FROMKEYS_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/types/dict.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <type_traits>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace dict
  {

    template <class Iterable, class V = types::none_type>
    types::dict<typename std::remove_reference_t<Iterable>::value_type, V>
    fromkeys(Iterable &&iter, V const &v = builtins::None);

    DEFINE_FUNCTOR(pythonic::builtins::dict, fromkeys);
  } // namespace dict
} // namespace builtins
PYTHONIC_NS_END

#endif
