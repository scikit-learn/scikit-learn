#ifndef PYTHONIC_NUMPY_REMAINDER_HPP
#define PYTHONIC_NUMPY_REMAINDER_HPP

#include "pythonic/include/numpy/remainder.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/numpy_broadcast.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME remainder
#define NUMPY_NARY_FUNC_SYM wrapper::remainder
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
