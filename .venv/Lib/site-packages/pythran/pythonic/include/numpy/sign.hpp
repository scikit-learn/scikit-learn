#ifndef PYTHONIC_INCLUDE_NUMPY_SIGN_HPP
#define PYTHONIC_INCLUDE_NUMPY_SIGN_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

#include <xsimd/xsimd.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {
#if PyArray_RUNTIME_VERSION < NPY_2_0_API_VERSION
    template <typename T>
    std::complex<T> sign(std::complex<T> v)
    {
      return xsimd::select(v == 0, v, v / xsimd::abs(v));
    }
#endif
    template <typename T>
    T sign(T v)
    {
      return xsimd::sign(v);
    }

  } // namespace details
#define NUMPY_NARY_FUNC_NAME sign
#define NUMPY_NARY_FUNC_SYM details::sign
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
