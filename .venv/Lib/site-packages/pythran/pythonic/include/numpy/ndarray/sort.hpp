#ifndef PYTHONIC_INCLUDE_NUMPY_NDARRAY_SORT_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDARRAY_SORT_HPP

#include "pythonic/include/numpy/sort.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace ndarray
  {
    template <class E>
    types::none_type sort(E &&expr, types::none_type);

    template <class E>
    types::none_type sort(E &&expr, long axis, types::none_type = {});

    template <class E>
    types::none_type sort(E &&expr, long axis, types::str const &kind);

    DEFINE_FUNCTOR(pythonic::numpy::ndarray, sort);
  } // namespace ndarray
} // namespace numpy
PYTHONIC_NS_END

#endif
