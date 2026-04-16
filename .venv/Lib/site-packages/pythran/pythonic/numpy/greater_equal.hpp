#ifndef PYTHONIC_NUMPY_GREATEREQUAL_HPP
#define PYTHONIC_NUMPY_GREATEREQUAL_HPP

#include "pythonic/include/numpy/greater_equal.hpp"

#include "pythonic/operator_/ge.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/numpy_broadcast.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME greater_equal
#define NUMPY_NARY_FUNC_SYM pythonic::operator_::ge
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
