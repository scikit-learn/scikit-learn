#ifndef PYTHONIC_INCLUDE_NUMPY_NANSUM_HPP
#define PYTHONIC_INCLUDE_NUMPY_NANSUM_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class F>
  void _nansum(E begin, E end, F &sum, utils::int_<1>);

  template <class E, class F, size_t N>
  void _nansum(E begin, E end, F &sum, utils::int_<N>);

  template <class E>
  typename E::dtype nansum(E const &expr);

  DEFINE_FUNCTOR(pythonic::numpy, nansum);
} // namespace numpy
PYTHONIC_NS_END

#endif
