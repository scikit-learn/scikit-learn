#ifndef PYTHONIC_NUMPY_SQRT_HPP
#define PYTHONIC_NUMPY_SQRT_HPP

#include "pythonic/include/numpy/sqrt.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME sqrt
#define NUMPY_NARY_FUNC_SYM xsimd::sqrt
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
