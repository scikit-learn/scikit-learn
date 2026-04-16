#ifndef PYTHONIC_INCLUDE_NUMPY_ARCTAN2_HPP
#define PYTHONIC_INCLUDE_NUMPY_ARCTAN2_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/numpy_broadcast.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

#include <xsimd/xsimd.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME arctan2
#define NUMPY_NARY_FUNC_SYM xsimd::atan2
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
