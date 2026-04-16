#ifndef PYTHONIC_NUMPY_ARCCOSH_HPP
#define PYTHONIC_NUMPY_ARCCOSH_HPP

#include "pythonic/include/numpy/arccosh.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME arccosh
#define NUMPY_NARY_FUNC_SYM xsimd::acosh
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
