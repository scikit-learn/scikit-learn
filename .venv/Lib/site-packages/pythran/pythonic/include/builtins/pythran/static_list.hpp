#ifndef PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_STATIC_LIST_HPP
#define PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_STATIC_LIST_HPP

#include "pythonic/include/builtins/list.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {
    inline types::empty_list static_list(std::tuple<> const &other)
    {
      return {};
    }
    template <class T, size_t N>
    types::static_list<T, N> static_list(types::array_tuple<T, N> const &other);
    template <class T, size_t N>
    types::static_list<T, N> static_list(types::array_tuple<T, N> &other);
    template <class T, size_t N>
    types::static_list<T, N> static_list(types::array_tuple<T, N> &&other);

    template <class T>
    auto static_list(T &&other)
        -> decltype(pythonic::builtins::functor::list{}(std::forward<T>(other)));

    template <class T0, class... Tys>
    types::static_list<typename __combined<T0, Tys...>::type, 1 + sizeof...(Tys)>
    static_list(std::tuple<T0, Tys...> const &other)
    {
      return static_list(types::to_array<typename __combined<T0, Tys...>::type>(other));
    }

    DEFINE_FUNCTOR(pythonic::builtins::pythran, static_list);
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
