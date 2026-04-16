#ifndef PYTHONIC_INCLUDE_BUILTIN_MINMAX_HPP
#define PYTHONIC_INCLUDE_BUILTIN_MINMAX_HPP

#include "pythonic/include/builtins/pythran/kwonly.hpp"
#include <utility>

PYTHONIC_NS_BEGIN

namespace builtins
{
  namespace details
  {
    template <class Op, class T>
    typename std::decay_t<T>::value_type minmax(Op const &, T &&t);

    template <class Op, class T, class F>
    typename std::decay_t<T>::value_type minmax(Op const &, T &&t, types::kwonly, F key);

    template <class Op, class T0, class T1, class... Types>
    std::enable_if_t<!std::is_same<T1, types::kwonly>::value,
                     typename __combined<T0, T1, Types...>::type>
    minmax(Op const &, T0 const &, T1 const &, Types const &...);
  } // namespace details
} // namespace builtins
PYTHONIC_NS_END

#endif
