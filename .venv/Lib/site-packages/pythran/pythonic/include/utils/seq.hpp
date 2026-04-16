#ifndef PYTHONIC_INCLUDE_UTILS_SEQ_HPP
#define PYTHONIC_INCLUDE_UTILS_SEQ_HPP

#include <utility>

PYTHONIC_NS_BEGIN

namespace utils
{

  namespace details
  {

    template <class T, T... Is>
    constexpr std::integer_sequence<T, (sizeof...(Is) - 1 - Is)...>
        reverse_integer_sequence(std::integer_sequence<T, Is...>);

  } // namespace details

  template <class T, std::size_t N>
  using make_reversed_integer_sequence =
      decltype(details::reverse_integer_sequence(std::make_integer_sequence<T, N>()));
  template <std::size_t N>
  using make_reversed_index_sequence = make_reversed_integer_sequence<std::size_t, N>;

  // make_repeated_type<A, 3>() => type_sequence<A, A, A>
  template <class... Tys>
  struct type_sequence {
  };

  namespace details
  {
    template <class T, std::size_t N, class... Tys>
    struct repeated_type : repeated_type<T, N - 1, T, Tys...> {
    };

    template <class T, class... Tys>
    struct repeated_type<T, 0, Tys...> {
      using type = type_sequence<Tys...>;
    };
  } // namespace details
  template <class T, std::size_t N>
  struct repeated_type : details::repeated_type<T, N> {
  };

  template <class T, std::size_t N>
  using make_repeated_type = typename repeated_type<T, N>::type;
} // namespace utils
PYTHONIC_NS_END

#endif
