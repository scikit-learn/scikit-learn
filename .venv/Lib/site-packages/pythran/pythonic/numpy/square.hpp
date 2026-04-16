#ifndef PYTHONIC_NUMPY_SQUARE_HPP
#define PYTHONIC_NUMPY_SQUARE_HPP

#include "pythonic/include/numpy/square.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME square
#define NUMPY_NARY_FUNC_SYM wrapper::square
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
