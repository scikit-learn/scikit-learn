#ifndef PYTHONIC_NUMPY_RECIPROCAL_HPP
#define PYTHONIC_NUMPY_RECIPROCAL_HPP

#include "pythonic/include/numpy/reciprocal.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME reciprocal
#define NUMPY_NARY_FUNC_SYM wrapper::reciprocal
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
