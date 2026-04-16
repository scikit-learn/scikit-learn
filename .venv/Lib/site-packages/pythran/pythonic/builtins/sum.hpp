#ifndef PYTHONIC_BUILTIN_SUM_HPP
#define PYTHONIC_BUILTIN_SUM_HPP

#include "pythonic/include/builtins/sum.hpp"

#include "pythonic/types/assignable.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/int_.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace details
  {
    template <class Tuple, size_t N>
    auto tuple_sum<Tuple, N>::operator()(Tuple const &t)
        -> decltype(std::get<N>(t) + tuple_sum<Tuple, N - 1>()(t))
    {
      return std::get<N>(t) + tuple_sum<Tuple, N - 1>()(t);
    }

    template <class Tuple>
    auto tuple_sum<Tuple, 0>::operator()(Tuple const &t) -> decltype(std::get<0>(t))
    {
      return std::get<0>(t);
    }
  } // namespace details

  template <class Iterable, class T>
  auto sum(Iterable s, T start) -> decltype(std::accumulate(
      s.begin(), s.end(),
      static_cast<typename assignable<decltype(start + *s.begin())>::type>(start)))
  {
    return std::accumulate(
        s.begin(), s.end(),
        static_cast<typename assignable<decltype(start + *s.begin())>::type>(start));
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
