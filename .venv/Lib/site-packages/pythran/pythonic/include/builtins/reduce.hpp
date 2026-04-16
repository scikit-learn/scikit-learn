#ifndef PYTHONIC_INCLUDE_BUILTIN_REDUCE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_REDUCE_HPP

#include "pythonic/include/utils/functor.hpp"

#include <numeric>
#include <utility>

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class Iterable, class Operator>
  auto reduce(Operator op, Iterable s) -> decltype(op(
      std::declval<typename std::iterator_traits<typename Iterable::iterator>::value_type>(),
      std::declval<typename std::iterator_traits<typename Iterable::iterator>::value_type>()));

  // this convoluted expression computes the fixed-point type of the output
  // it's required because, e.g. static_list<long, 1> + static_list<long, 1>
  // returns array_tuple<long, 2>
  // and this widens to list
  template <class Iterable, class Operator, class T>
  using reduce_helper_t =
      typename __combined<T, decltype(std::declval<Operator>()(
                                 std::declval<T const &>(),
                                 std::declval<typename std::iterator_traits<
                                     typename Iterable::iterator>::value_type>()))>::type;

  template <class Iterable, class Operator, class T>
  auto reduce(Operator op, Iterable s, T const &init)
      -> decltype(std::accumulate(s.begin(), s.end(),
                                  static_cast<reduce_helper_t<Iterable, Operator, T>>(init), op));

  DEFINE_FUNCTOR(pythonic::builtins, reduce);
} // namespace builtins
PYTHONIC_NS_END

#endif
