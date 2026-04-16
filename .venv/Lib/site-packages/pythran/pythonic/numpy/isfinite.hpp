#ifndef PYTHONIC_NUMPY_ISFINITE_HPP
#define PYTHONIC_NUMPY_ISFINITE_HPP

#include "pythonic/include/numpy/isfinite.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME isfinite
#define NUMPY_NARY_FUNC_SYM wrapper::isfinite
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
