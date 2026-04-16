#ifndef PYTHONIC_INCLUDE_NUMPY_PLACE_HPP
#define PYTHONIC_INCLUDE_NUMPY_PLACE_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS, class Tp, class pSp, class F>
  types::none_type place(types::ndarray<T, pS> &expr, types::ndarray<Tp, pSp> const &mask,
                         F const &values);

  template <class T, class pS, class M, class F>
  types::none_type place(types::ndarray<T, pS> &expr, M const &mask, F const &values);

  template <class E, class M, class F>
  types::none_type place(E &, M const &, F const &);

  DEFINE_FUNCTOR(pythonic::numpy, place);
} // namespace numpy
PYTHONIC_NS_END

#endif
