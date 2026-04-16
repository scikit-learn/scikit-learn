#ifndef PYTHONIC_INCLUDE_NUMPY_BITWISE_AND_HPP
#define PYTHONIC_INCLUDE_NUMPY_BITWISE_AND_HPP

#include "pythonic/include/operator_/and_.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/numpy_broadcast.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME bitwise_and
#define NUMPY_NARY_FUNC_SYM pythonic::operator_::and_
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
