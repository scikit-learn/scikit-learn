#ifndef PYTHONIC_INCLUDE_NUMPY_LEFT_SHIFT_HPP
#define PYTHONIC_INCLUDE_NUMPY_LEFT_SHIFT_HPP

#include "pythonic/include/operator_/lshift.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/numpy_broadcast.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME left_shift
#define NUMPY_NARY_FUNC_SYM pythonic::operator_::lshift
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
