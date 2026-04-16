#ifndef PYTHONIC_INCLUDE_NUMPY_RESIZE_HPP
#define PYTHONIC_INCLUDE_NUMPY_RESIZE_HPP

#include "pythonic/include/numpy/ndarray/reshape.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  USING_FUNCTOR(resize, pythonic::numpy::ndarray::functor::reshape);
}
PYTHONIC_NS_END

#endif
