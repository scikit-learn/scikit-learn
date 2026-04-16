#ifndef PYTHONIC_NUMPY_FIX_HPP
#define PYTHONIC_NUMPY_FIX_HPP

#include "pythonic/include/numpy/fix.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME fix
#define NUMPY_NARY_FUNC_SYM wrapper::fix
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
