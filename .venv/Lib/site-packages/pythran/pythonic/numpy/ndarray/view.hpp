#ifndef PYTHONIC_NUMPY_NDARRAY_VIEW_HPP
#define PYTHONIC_NUMPY_NDARRAY_VIEW_HPP

#include "pythonic/include/numpy/ndarray/view.hpp"
#include "pythonic/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace ndarray
  {
    template <class E>
    auto view(E &&e) -> decltype(std::forward<E>(e))
    {
      return std::forward<E>(e);
    }

    template <class E, class dtype>
    auto view(E &&e, dtype d)
        -> decltype(std::forward<E>(e).template recast<typename dtype::type>())
    {
      return std::forward<E>(e).template recast<typename dtype::type>();
    }

  } // namespace ndarray
} // namespace numpy
PYTHONIC_NS_END

#endif
