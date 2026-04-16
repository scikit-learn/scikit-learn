#ifndef PYTHONIC_INCLUDE_NUMPY_UNION1D_HPP
#define PYTHONIC_INCLUDE_NUMPY_UNION1D_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class F>
  types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                 types::pshape<long>>
  union1d(E const &e, F const &f);

  DEFINE_FUNCTOR(pythonic::numpy, union1d)
} // namespace numpy
PYTHONIC_NS_END

#endif
