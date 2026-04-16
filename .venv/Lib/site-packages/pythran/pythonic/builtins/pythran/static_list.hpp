#ifndef PYTHONIC_BUILTIN_PYTHRAN_STATIC_LIST_HPP
#define PYTHONIC_BUILTIN_PYTHRAN_STATIC_LIST_HPP

#include "pythonic/builtins/list.hpp"
#include "pythonic/include/builtins/pythran/static_list.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {
    template <class T, size_t N>
    types::static_list<T, N> static_list(types::array_tuple<T, N> const &other)
    {
      return other.template to_array<types::list_version>();
    }
    template <class T, size_t N>
    types::static_list<T, N> static_list(types::array_tuple<T, N> &other)
    {
      return other.template to_array<types::list_version>();
    }
    template <class T, size_t N>
    types::static_list<T, N> static_list(types::array_tuple<T, N> &&other)
    {
      return other.template to_array<types::list_version>();
    }

    template <class T>
    auto static_list(T &&other)
        -> decltype(pythonic::builtins::functor::list{}(std::forward<T>(other)))
    {
      return pythonic::builtins::functor::list{}(std::forward<T>(other));
    }
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
