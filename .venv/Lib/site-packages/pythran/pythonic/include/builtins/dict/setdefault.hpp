#ifndef PYTHONIC_INCLUDE_BUILTIN_DICT_SETDEFAULT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_DICT_SETDEFAULT_HPP

#include "pythonic/include/types/dict.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace dict
  {

    template <class K, class V, class W, class X>
    V &setdefault(types::dict<K, V> &d, W const &k, X const &default_);

    template <class K, class V, class W>
    types::none<V> setdefault(types::dict<K, V> &d, W const &k);

    template <class K, class V, class W, class X>
    V setdefault(types::dict<K, V> &&d, W const &k, X const &default_);

    template <class K, class V, class W>
    types::none<V> setdefault(types::dict<K, V> &&d, W const &k);

    DEFINE_FUNCTOR(pythonic::builtins::dict, setdefault);
  } // namespace dict
} // namespace builtins
PYTHONIC_NS_END

#endif
