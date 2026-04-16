#ifndef PYTHONIC_NUMPY_NEXTAFTER_HPP
#define PYTHONIC_NUMPY_NEXTAFTER_HPP

#include "pythonic/include/numpy/nextafter.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/numpy_broadcast.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
#define NUMPY_NARY_FUNC_NAME nextafter
#define NUMPY_NARY_FUNC_SYM std::nextafter
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
