#ifndef PYTHONIC_INCLUDE_NUMPY_BROADCAST_TO_HPP
#define PYTHONIC_INCLUDE_NUMPY_BROADCAST_TO_HPP

#include "pythonic/include/numpy/empty.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class pS>
  auto broadcast_to(E const &expr, pS shape) -> decltype(numpy::functor::empty{}(
      shape, typename types::dtype_t<typename types::dtype_of<E>::type>{}));

  DEFINE_FUNCTOR(pythonic::numpy, broadcast_to);
} // namespace numpy
PYTHONIC_NS_END

#endif
