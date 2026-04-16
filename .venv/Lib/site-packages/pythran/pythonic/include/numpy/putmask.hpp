#ifndef PYTHONIC_INCLUDE_NUMPY_PUTMASK_HPP
#define PYTHONIC_INCLUDE_NUMPY_PUTMASK_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS, class E, class F>
  types::none_type putmask(types::ndarray<T, pS> &expr, E const &mask, F const &values);

  template <class E, class M, class F>
  types::none_type putmask(E &, M const &, F const &);

  DEFINE_FUNCTOR(pythonic::numpy, putmask);
} // namespace numpy
PYTHONIC_NS_END

#endif
