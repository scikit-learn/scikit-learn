#ifndef PYTHONIC_NUMPY_RIGHTSHIFT_HPP
#define PYTHONIC_NUMPY_RIGHTSHIFT_HPP

#include "pythonic/include/numpy/right_shift.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/numpy_broadcast.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

#include "pythonic/operator_/rshift.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME right_shift
#define NUMPY_NARY_FUNC_SYM operator_::rshift
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
