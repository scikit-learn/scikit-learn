#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_MOD_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_MOD_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace str
  {
    template <class T>
    types::str __mod__(types::str const &, T const &arg);
    template <class... Ts>
    types::str __mod__(types::str const &, std::tuple<Ts...> const &args);
    template <size_t N, class T>
    types::str __mod__(types::str const &, types::array_tuple<T, N> const &args);

    DEFINE_FUNCTOR(pythonic::builtins::str, __mod__);
  } // namespace str
} // namespace builtins
PYTHONIC_NS_END
#endif
