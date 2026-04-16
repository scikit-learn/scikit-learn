#ifndef PYTHONIC_NUMPY_PUTMASK_HPP
#define PYTHONIC_NUMPY_PUTMASK_HPP

#include "pythonic/include/numpy/putmask.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS, class E, class F>
  types::none_type putmask(types::ndarray<T, pS> &expr, E const &mask, F const &values)
  {
    auto amask = asarray(mask);
    auto avalues = asarray(values);
    auto iexpr = expr.fbegin();
    auto n = avalues.flat_size();
    for (long i = 0; i < expr.flat_size(); ++i)
      if (*(amask.fbegin() + i))
        *(iexpr + i) = *(avalues.fbegin() + i % n);
    return builtins::None;
  }

  template <class E, class M, class F>
  types::none_type putmask(E &, M const &, F const &)
  {
    throw std::runtime_error("putmask only partially implemented");
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
