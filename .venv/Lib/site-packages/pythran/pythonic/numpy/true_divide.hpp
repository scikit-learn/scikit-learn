#ifndef PYTHONIC_NUMPY_TRUEDIVIDE_HPP
#define PYTHONIC_NUMPY_TRUEDIVIDE_HPP

#include "pythonic/include/numpy/true_divide.hpp"

#include "pythonic/operator_/div.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/numpy_broadcast.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

// FIXME: this is ! always a true_divide...
#define NUMPY_NARY_FUNC_NAME true_divide
#define NUMPY_NARY_FUNC_SYM pythonic::operator_::div
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
