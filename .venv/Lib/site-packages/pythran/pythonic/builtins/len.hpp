#ifndef PYTHONIC_BUILTIN_LEN_HPP
#define PYTHONIC_BUILTIN_LEN_HPP

#include "pythonic/include/builtins/len.hpp"

#include "pythonic/types/traits.hpp"
#include "pythonic/utils/functor.hpp"

#include <iterator>
#include <tuple>

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class... Types>
  long len(std::tuple<Types...> const &)
  {
    return sizeof...(Types);
  }

  template <class T>
  std::enable_if_t<types::has_size<T>::value, long> len(T const &t)
  {
    return t.size();
  }
} // namespace builtins
PYTHONIC_NS_END
#endif
