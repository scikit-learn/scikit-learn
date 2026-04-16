#ifndef PYTHONIC_INCLUDE_BUILTIN_SUM_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SUM_HPP

#include "pythonic/include/types/assignable.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/int_.hpp"

#include <numeric>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace details
  {
    template <class Tuple, size_t N>
    struct tuple_sum {
      auto operator()(Tuple const &t) -> decltype(std::get<N>(t) + tuple_sum<Tuple, N - 1>()(t));
    };

    template <class Tuple>
    struct tuple_sum<Tuple, 0> {
      auto operator()(Tuple const &t) -> decltype(std::get<0>(t));
    };
  } // namespace details

  template <class Iterable, class T>
  auto sum(Iterable s, T start) -> decltype(std::accumulate(
      s.begin(), s.end(),
      static_cast<typename assignable<decltype(start + *s.begin())>::type>(start)));

  template <class Iterable>
  auto sum(Iterable s) -> decltype(sum(s, 0L))
  {
    return sum(s, 0L);
  }

  template <class... Types>
  auto sum(std::tuple<Types...> const &t)
      -> decltype(details::tuple_sum<std::tuple<Types...>, sizeof...(Types) - 1>()(t))
  {
    return details::tuple_sum<std::tuple<Types...>, sizeof...(Types) - 1>()(t);
  }

  template <class T, size_t N, class V>
  T sum(types::array_base<T, N, V> const &t)
  {
    return details::tuple_sum<types::array_base<T, N, V>, N - 1>()(t);
  }

  DEFINE_FUNCTOR(pythonic::builtins, sum);
} // namespace builtins
PYTHONIC_NS_END

#endif
