#ifndef PYTHONIC_NUMPY_COPYSIGN_HPP
#define PYTHONIC_NUMPY_COPYSIGN_HPP

#include "pythonic/include/numpy/copysign.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/numpy_broadcast.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME copysign
#define NUMPY_NARY_FUNC_SYM xsimd::copysign
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
