#ifndef PYTHONIC_INCLUDE_NUMPY_FLOOR_HPP
#define PYTHONIC_INCLUDE_NUMPY_FLOOR_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

#include <xsimd/xsimd.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME floor
#define NUMPY_NARY_FUNC_SYM xsimd::floor
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
