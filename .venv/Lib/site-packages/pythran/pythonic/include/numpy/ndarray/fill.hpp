#ifndef PYTHONIC_INCLUDE_NUMPY_NDARRAY_FILL_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDARRAY_FILL_HPP

#include "pythonic/builtins/None.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace ndarray
  {
    template <class E, class F>
    types::none_type fill(E &&e, F f);

    template <class T, class pS, class F>
    types::none_type fill(types::ndarray<T, pS> &e, F f);

    DEFINE_FUNCTOR(pythonic::numpy::ndarray, fill);
  } // namespace ndarray
} // namespace numpy
PYTHONIC_NS_END

#endif
