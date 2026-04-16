#ifndef PYTHONIC_NUMPY_ASARRAYCHKFINITE_HPP
#define PYTHONIC_NUMPY_ASARRAYCHKFINITE_HPP

#include "pythonic/include/numpy/asarray_chkfinite.hpp"

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/numpy/isfinite.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace wrapper
  {
    template <class I>
    I asarray_chkfinite(I const &a)
    {
      if (!functor::isfinite()(a))
        throw types::ValueError("array must ! contain infs || NaNs");
      return a;
    }
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME asarray_chkfinite
#define NUMPY_NARY_FUNC_SYM wrapper::asarray_chkfinite
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
