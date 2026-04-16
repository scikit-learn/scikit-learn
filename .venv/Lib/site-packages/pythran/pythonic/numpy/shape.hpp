#ifndef PYTHONIC_NUMPY_SHAPE_HPP
#define PYTHONIC_NUMPY_SHAPE_HPP

#include "pythonic/include/numpy/shape.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class T, class pS>
  auto shape(types::ndarray<T, pS> const &e) -> decltype(e._shape)
  {
    return e._shape;
  }

  template <class E>
  auto shape(E const &e) -> decltype(sutils::getshape(e))
  {
    return sutils::getshape(e);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
