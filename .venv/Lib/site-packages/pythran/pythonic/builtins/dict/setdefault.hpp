#ifndef PYTHONIC_BUILTIN_DICT_SETDEFAULT_HPP
#define PYTHONIC_BUILTIN_DICT_SETDEFAULT_HPP

#include "pythonic/include/builtins/dict/setdefault.hpp"

#include "pythonic/types/dict.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace dict
  {

    template <class K, class V, class W, class X>
    V &setdefault(types::dict<K, V> &d, W const &k, X const &default_)
    {
      return d.setdefault(k, default_);
    }

    template <class K, class V, class W>
    types::none<V> setdefault(types::dict<K, V> &d, W const &k)
    {
      return d.get(k);
    }

    template <class K, class V, class W, class X>
    V setdefault(types::dict<K, V> &&d, W const &k, X const &default_)
    {
      return d.setdefault(k, default_);
    }

    template <class K, class V, class W>
    types::none<V> setdefault(types::dict<K, V> &&d, W const &k)
    {
      return d.get(k);
    }
  } // namespace dict
} // namespace builtins
PYTHONIC_NS_END

#endif
