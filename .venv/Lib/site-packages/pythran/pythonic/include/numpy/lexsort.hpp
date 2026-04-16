#ifndef PYTHONIC_INCLUDE_NUMPY_LEXSORT_HPP
#define PYTHONIC_INCLUDE_NUMPY_LEXSORT_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class pS>
  types::ndarray<long, types::pshape<long>> lexsort(pS const &keys);

  DEFINE_FUNCTOR(pythonic::numpy, lexsort)
} // namespace numpy
PYTHONIC_NS_END

#endif
