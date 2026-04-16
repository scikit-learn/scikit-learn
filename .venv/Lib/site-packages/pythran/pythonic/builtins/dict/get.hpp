#ifndef PYTHONIC_BUILTIN_DICT_GET_HPP
#define PYTHONIC_BUILTIN_DICT_GET_HPP

#include "pythonic/include/builtins/dict/get.hpp"

#include "pythonic/types/NoneType.hpp"
#include "pythonic/types/dict.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace dict
  {
    template <class K, class V, class W, class X>
    typename __combined<V, X>::type get(types::dict<K, V> const &d, W const &k, X const &default_)
    {
      return d.get(k, default_);
    }

    template <class K, class V, class W>
    types::none<V> get(types::dict<K, V> const &d, W const &k)
    {
      return d.get(k);
    }

    template <class W, class X>
    X get(types::empty_dict const &, W const &, X const &default_)
    {
      return default_;
    }
  } // namespace dict
} // namespace builtins
PYTHONIC_NS_END

#endif
