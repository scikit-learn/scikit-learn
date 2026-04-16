#ifndef PYTHONIC_NUMPY_NEGATIVE_HPP
#define PYTHONIC_NUMPY_NEGATIVE_HPP

#include "pythonic/include/numpy/negative.hpp"

#include "pythonic/operator_/neg.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME negative
#define NUMPY_NARY_FUNC_SYM pythonic::operator_::neg
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
