#ifndef PYTHONIC_INCLUDE_NUMPY_RINT_HPP
#define PYTHONIC_INCLUDE_NUMPY_RINT_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    T rint(T const &v);
    template <class T>
    std::complex<T> rint(std::complex<T> const &v);
  } // namespace wrapper
#define NUMPY_NARY_FUNC_NAME rint
#define NUMPY_NARY_FUNC_SYM wrapper::rint
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
