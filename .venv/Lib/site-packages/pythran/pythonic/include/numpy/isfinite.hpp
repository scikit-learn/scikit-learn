#ifndef PYTHONIC_INCLUDE_NUMPY_ISFINITE_HPP
#define PYTHONIC_INCLUDE_NUMPY_ISFINITE_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    bool isfinite(std::complex<T> const &t)
    {
      return std::isfinite(t.real()) && std::isfinite(t.imag());
    }
    template <class T>
    bool isfinite(T const &v)
    {
      return std::isfinite(v);
    }
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME isfinite
#define NUMPY_NARY_FUNC_SYM wrapper::isfinite
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
