#ifndef PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_ABSSQR_HPP
#define PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_ABSSQR_HPP

#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace pythran
  {
    namespace details
    {

      template <class T>
      T abssqr(T const &v);

      template <class T>
      T abssqr(std::complex<T> const &v);
    } // namespace details

#define NUMPY_NARY_FUNC_NAME abssqr
#define NUMPY_NARY_FUNC_SYM details::abssqr
#include "pythonic/include/types/numpy_nary_expr.hpp"
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
