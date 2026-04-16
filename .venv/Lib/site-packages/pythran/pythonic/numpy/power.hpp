#ifndef PYTHONIC_NUMPY_POWER_HPP
#define PYTHONIC_NUMPY_POWER_HPP

#include "pythonic/include/numpy/power.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME power
#define NUMPY_NARY_FUNC_SYM xsimd::pow
// no need to adapt_type here, as it may turn a**2 into a**2.f
#define NUMPY_NARY_RESHAPE_MODE reshape_type
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
