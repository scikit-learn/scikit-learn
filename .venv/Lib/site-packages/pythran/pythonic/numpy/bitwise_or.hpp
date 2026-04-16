#ifndef PYTHONIC_NUMPY_BITWISE_OR_HPP
#define PYTHONIC_NUMPY_BITWISE_OR_HPP

#include "pythonic/include/numpy/bitwise_or.hpp"

#include "pythonic/operator_/or_.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/numpy_broadcast.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME bitwise_or
#define NUMPY_NARY_FUNC_SYM pythonic::operator_::or_
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
