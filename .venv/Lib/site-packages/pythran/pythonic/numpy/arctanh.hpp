#ifndef PYTHONIC_NUMPY_ARCTANH_HPP
#define PYTHONIC_NUMPY_ARCTANH_HPP

#include "pythonic/include/numpy/arctanh.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME arctanh
#define NUMPY_NARY_FUNC_SYM xsimd::atanh
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
