#ifndef PYTHONIC_NUMPY_RAD2DEG_HPP
#define PYTHONIC_NUMPY_RAD2DEG_HPP

#include "pythonic/include/numpy/rad2deg.hpp"

#include "pythonic/numpy/pi.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME rad2deg
#define NUMPY_NARY_FUNC_SYM wrapper::rad2deg
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
