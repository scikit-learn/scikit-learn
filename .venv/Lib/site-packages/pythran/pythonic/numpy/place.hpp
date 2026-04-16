#ifndef PYTHONIC_NUMPY_PLACE_HPP
#define PYTHONIC_NUMPY_PLACE_HPP

#include "pythonic/include/numpy/place.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS, class Tp, class pSp, class F>
  types::none_type place(types::ndarray<T, pS> &expr, types::ndarray<Tp, pSp> const &mask,
                         F const &values)
  {
    auto first = expr.fend();
    auto viter = values.begin(), vend = values.end();
    auto miter = mask.fbegin();
    for (auto iter = expr.fbegin(), end = expr.fend(); iter != end; ++iter, ++miter) {
      if (*miter) {
        if (first == expr.fend())
          first = iter;
        if (viter == vend)
          viter = values.begin();
        *iter = *viter;
        ++viter;
      }
    }
    return builtins::None;
  }

  template <class T, class pS, class M, class F>
  types::none_type place(types::ndarray<T, pS> &expr, M const &mask, F const &values)
  {
    return place(expr, asarray(mask), values);
  }

  template <class E, class M, class F>
  types::none_type place(E &, M const &, F const &)
  {
    throw std::runtime_error("place only partially implemented");
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
