#ifndef PYTHONIC_NUMPY_LESS_HPP
#define PYTHONIC_NUMPY_LESS_HPP

#include "pythonic/include/numpy/less.hpp"

#include "pythonic/operator_/lt.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/numpy_broadcast.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME less
#define NUMPY_NARY_FUNC_SYM pythonic::operator_::lt
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
