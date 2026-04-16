#ifndef PYTHONIC_INCLUDE_NUMPY_DOUBLE_HPP
#define PYTHONIC_INCLUDE_NUMPY_DOUBLE_HPP

#include "pythonic/include/numpy/float64.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME double_
#define NUMPY_NARY_FUNC_SYM details::float64
#define NUMPY_NARY_EXTRA_METHOD using type = double;
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
