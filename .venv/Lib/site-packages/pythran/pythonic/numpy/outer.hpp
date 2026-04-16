#ifndef PYTHONIC_NUMPY_OUTER_HPP
#define PYTHONIC_NUMPY_OUTER_HPP

#include "pythonic/include/numpy/outer.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T0, class pS0, class T1, class pS1>
  types::ndarray<decltype(std::declval<T0>() + std::declval<T1>()), types::pshape<long, long>>
  outer(types::ndarray<T0, pS0> const &a, types::ndarray<T1, pS1> const &b)
  {
    types::ndarray<decltype(std::declval<T0>() + std::declval<T1>()), types::pshape<long, long>>
        out(types::pshape<long, long>{a.flat_size(), b.flat_size()}, builtins::None);
    auto iter = out.fbegin();
    for (auto iter_a = a.fbegin(), end_a = a.fend(); iter_a != end_a; ++iter_a) {
      auto val_a = *iter_a;
      iter = std::transform(b.fbegin(), b.fend(), iter, [=](T1 val) { return val_a * val; });
    }
    return out;
  }

  template <class T0, class pS0, class E1>
  auto outer(types::ndarray<T0, pS0> const &a, E1 const &b) -> decltype(outer(a, asarray(b)))
  {
    return outer(a, asarray(b));
  }

  template <class E0, class T1, class pS1>
  auto outer(E0 const &a, types::ndarray<T1, pS1> const &b) -> decltype(outer(asarray(a), b))
  {
    return outer(asarray(a), b);
  }

  template <class E0, class E1>
  auto outer(E0 const &a, E1 const &b) -> decltype(outer(asarray(a), asarray(b)))
  {
    return outer(asarray(a), asarray(b));
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
