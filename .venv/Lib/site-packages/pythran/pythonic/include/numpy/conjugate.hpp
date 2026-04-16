#ifndef PYTHONIC_INCLUDE_NUMPY_CONJUGATE_HPP
#define PYTHONIC_INCLUDE_NUMPY_CONJUGATE_HPP

#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

#include <xsimd/xsimd.hpp>

// Inject some extra symbol in xsimd namespace, until that's fixed upstream
#if XSIMD_VERSION_MAJOR < 8 ||                                                                     \
    (XSIMD_VERSION_MAJOR == 8 && XSIMD_VERSION_MINOR == 0 && XSIMD_VERSION_PATCH <= 5)
namespace xsimd
{
  using std::conj;
}
#endif

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    std::complex<T> conjugate(std::complex<T> const &v)
    {
      return std::conj(v);
    }
    template <class T, class A>
    xsimd::batch<std::complex<T>, A> conjugate(xsimd::batch<std::complex<T>, A> const &v)
    {
      return xsimd::conj(v);
    }
    template <class T, class A>
    xsimd::batch<T, A> conjugate(xsimd::batch<T, A> const &v)
    {
      return v;
    }
    template <class T>
    T conjugate(T const &v)
    {
      return v;
    }
  } // namespace wrapper
#define NUMPY_NARY_FUNC_NAME conjugate
#define NUMPY_NARY_FUNC_SYM wrapper::conjugate
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
