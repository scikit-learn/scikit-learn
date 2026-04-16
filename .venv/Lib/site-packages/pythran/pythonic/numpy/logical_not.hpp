#ifndef PYTHONIC_NUMPY_LOGICALNOT_HPP
#define PYTHONIC_NUMPY_LOGICALNOT_HPP

#include "pythonic/include/numpy/logical_not.hpp"

#include "pythonic/operator_/not_.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME logical_not
#define NUMPY_NARY_FUNC_SYM pythonic::operator_::not_
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
