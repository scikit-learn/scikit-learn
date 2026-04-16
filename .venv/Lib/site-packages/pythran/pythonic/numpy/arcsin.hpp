#ifndef PYTHONIC_NUMPY_ARCSIN_HPP
#define PYTHONIC_NUMPY_ARCSIN_HPP

#include "pythonic/include/numpy/arcsin.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME arcsin
#define NUMPY_NARY_FUNC_SYM xsimd::asin
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
