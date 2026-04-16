#ifndef PYTHONIC_NUMPY_RINT_HPP
#define PYTHONIC_NUMPY_RINT_HPP

#include "pythonic/include/numpy/rint.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    T rint(T const &v)
    {
      return std::nearbyint(v);
    }
    template <class T>
    std::complex<T> rint(std::complex<T> const &v)
    {
      return {std::nearbyint(v.real()), std::nearbyint(v.imag())};
    }
  } // namespace wrapper
#define NUMPY_NARY_FUNC_NAME rint
#define NUMPY_NARY_FUNC_SYM wrapper::rint
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
