#ifndef PYTHONIC_INCLUDE_NUMPY_NDARRAY_VIEW_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDARRAY_VIEW_HPP

#include "pythonic/include/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace ndarray
  {

    template <class E>
    auto view(E &&e) -> decltype(std::forward<E>(e));

    template <class E, class dtype>
    auto view(E &&e, dtype d)
        -> decltype(std::forward<E>(e).template recast<typename dtype::type>());

    DEFINE_FUNCTOR(pythonic::numpy::ndarray, view);
  } // namespace ndarray
} // namespace numpy
PYTHONIC_NS_END

#endif
