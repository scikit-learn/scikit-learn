#ifndef PYTHONIC_INCLUDE_NUMPY_SHAPE_HPP
#define PYTHONIC_INCLUDE_NUMPY_SHAPE_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class T, class pS>
  auto shape(types::ndarray<T, pS> const &e) -> decltype(e._shape);

  template <class E>
  auto shape(E const &e) -> decltype(sutils::getshape(e));

  DEFINE_FUNCTOR(pythonic::numpy, shape)
} // namespace numpy
PYTHONIC_NS_END

#endif
