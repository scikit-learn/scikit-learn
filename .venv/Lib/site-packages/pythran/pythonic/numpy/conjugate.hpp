#ifndef PYTHONIC_NUMPY_CONJUGATE_HPP
#define PYTHONIC_NUMPY_CONJUGATE_HPP

#include "pythonic/include/numpy/conjugate.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME conjugate
#define NUMPY_NARY_FUNC_SYM wrapper::conjugate
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
