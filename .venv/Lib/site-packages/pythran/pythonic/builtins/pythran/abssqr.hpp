#ifndef PYTHONIC_BUILTIN_PYTHRAN_ABSSQR_HPP
#define PYTHONIC_BUILTIN_PYTHRAN_ABSSQR_HPP

#include "pythonic/include/builtins/pythran/abssqr.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  namespace pythran
  {

    namespace details
    {

      template <class T>
      T abssqr(T const &v)
      {
        return v * v;
      }

      template <class T>
      T abssqr(std::complex<T> const &v)
      {
        return v.real() * v.real() + v.imag() * v.imag();
      }
    } // namespace details

#define NUMPY_NARY_FUNC_NAME abssqr
#define NUMPY_NARY_FUNC_SYM details::abssqr
#include "pythonic/types/numpy_nary_expr.hpp"
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
