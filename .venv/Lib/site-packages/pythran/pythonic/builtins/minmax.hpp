#ifndef PYTHONIC_BUILTIN_MINMAX_HPP
#define PYTHONIC_BUILTIN_MINMAX_HPP

#include "pythonic/include/builtins/minmax.hpp"

#include <algorithm>
#include <utility>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace details
  {
    template <typename Op, class T>
    typename std::decay_t<T>::value_type minmax(Op const &op, T &&t)
    {
      return *std::max_element(t.begin(), t.end(), op);
    }

    template <typename Op, class T, class F>
    typename std::decay_t<T>::value_type minmax(Op const &op, T &&t, types::kwonly, F key)
    {
      using value_type = decltype(*t.begin());
      return *std::max_element(t.begin(), t.end(),
                               [op, key](value_type const &self, value_type const &other) {
                                 return op(key(self), key(other));
                               });
    }

    template <class Op, class T0, class T1, class... Types>
    std::enable_if_t<!std::is_same<T1, types::kwonly>::value,
                     typename __combined<T0, T1, Types...>::type>
    minmax(Op const &op, T0 const &t0, T1 const &t1, Types const &...ts)
    {
      using value_type = typename __combined<T0, T1, Types...>::type;
      std::initializer_list<value_type> values = {
          static_cast<value_type>(t0), static_cast<value_type>(t1), static_cast<value_type>(ts)...};
      return minmax(op, values);
    }
  } // namespace details
} // namespace builtins
PYTHONIC_NS_END

#endif
