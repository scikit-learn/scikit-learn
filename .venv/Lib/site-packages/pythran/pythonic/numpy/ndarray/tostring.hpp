#ifndef PYTHONIC_NUMPY_NDARRAY_TOSTRING_HPP
#define PYTHONIC_NUMPY_NDARRAY_TOSTRING_HPP

#include "pythonic/include/numpy/ndarray/tostring.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace ndarray
  {
    template <class T, class pS>
    types::str tostring(types::ndarray<T, pS> const &expr)
    {
      return types::str(reinterpret_cast<const char *>(expr.buffer), expr.flat_size() * sizeof(T));
    }
    NUMPY_EXPR_TO_NDARRAY0_IMPL(tostring);
  } // namespace ndarray
} // namespace numpy
PYTHONIC_NS_END
#endif
