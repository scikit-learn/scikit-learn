#ifndef PYTHONIC_NUMPY_ANGLEINDEG_HPP
#define PYTHONIC_NUMPY_ANGLEINDEG_HPP

#include "pythonic/include/numpy/angle_in_deg.hpp"

#include "pythonic/numpy/angle_in_rad.hpp"
#include "pythonic/numpy/pi.hpp"
#include "pythonic/utils/numpy_traits.hpp"

/* NOTE: angle_in_deg is not part of the official Numpy API,
 * this file is here only to split the angle function in two parts
 */

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME angle_in_deg
#define NUMPY_NARY_FUNC_SYM wrapper::angle_in_deg
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
