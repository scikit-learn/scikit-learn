#ifndef PYTHONIC_INCLUDE_NUMPY_SETDIFF1D_HPP
#define PYTHONIC_INCLUDE_NUMPY_SETDIFF1D_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class T, class U>
  types::ndarray<typename __combined<typename types::dtype_of<T>::type,
                                     typename types::dtype_of<U>::type>::type,
                 types::pshape<long>>
  setdiff1d(T const &ar1, U const &ar2, bool assume_unique = false);

  DEFINE_FUNCTOR(pythonic::numpy, setdiff1d);
} // namespace numpy
PYTHONIC_NS_END

#endif
