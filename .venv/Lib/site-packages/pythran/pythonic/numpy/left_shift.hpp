#ifndef PYTHONIC_NUMPY_LEFT_SHIFT_HPP
#define PYTHONIC_NUMPY_LEFT_SHIFT_HPP

#include "pythonic/include/numpy/left_shift.hpp"

#include "pythonic/operator_/lshift.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/numpy_broadcast.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME left_shift
#define NUMPY_NARY_FUNC_SYM pythonic::operator_::lshift
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
