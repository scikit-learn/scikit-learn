#ifndef PYTHONIC_NUMPY_COSH_HPP
#define PYTHONIC_NUMPY_COSH_HPP

#include "pythonic/include/numpy/cosh.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

#include <xsimd/xsimd.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME cosh
#define NUMPY_NARY_FUNC_SYM xsimd::cosh
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
