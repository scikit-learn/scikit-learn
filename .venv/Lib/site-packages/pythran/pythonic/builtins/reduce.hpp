#ifndef PYTHONIC_BUILTIN_REDUCE_HPP
#define PYTHONIC_BUILTIN_REDUCE_HPP

#include "pythonic/include/builtins/reduce.hpp"

#include "pythonic/utils/functor.hpp"

#include <algorithm>
#include <numeric>
#include <utility>

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class Iterable, class Operator>
  auto reduce(Operator op, Iterable s) -> decltype(op(
      std::declval<typename std::iterator_traits<typename Iterable::iterator>::value_type>(),
      std::declval<typename std::iterator_traits<typename Iterable::iterator>::value_type>()))
  {
    auto iter = s.begin();
    auto r = *iter;
    ++iter;
    if (iter != s.end())
      return std::accumulate(iter, s.end(), r, op);
    else
      return r;
  }

  template <class Iterable, class Operator, class T>
  auto reduce(Operator op, Iterable s, T const &init)
      -> decltype(std::accumulate(s.begin(), s.end(),
                                  static_cast<reduce_helper_t<Iterable, Operator, T>>(init), op))
  {
    return std::accumulate(s.begin(), s.end(),
                           static_cast<reduce_helper_t<Iterable, Operator, T>>(init), op);
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
