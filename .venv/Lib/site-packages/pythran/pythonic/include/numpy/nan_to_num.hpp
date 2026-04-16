#ifndef PYTHONIC_INCLUDE_NUMPY_NANTONUM_HPP
#define PYTHONIC_INCLUDE_NUMPY_NANTONUM_HPP

#include "pythonic/include/numpy/isnan.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

#include <limits>

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace wrapper
  {
    template <class I>
    I nan_to_num(I const &a);
  }

#define NUMPY_NARY_FUNC_NAME nan_to_num
#define NUMPY_NARY_FUNC_SYM wrapper::nan_to_num
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
