#ifndef PYTHONIC_INCLUDE_NUMPY_SPACING_HPP
#define PYTHONIC_INCLUDE_NUMPY_SPACING_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    auto spacing(T const &v) -> decltype(std::nextafter(v, 1) - v)
    {
      return std::nextafter(v, 1) - v;
    }
  } // namespace wrapper
#define NUMPY_NARY_FUNC_NAME spacing
#define NUMPY_NARY_FUNC_SYM wrapper::spacing
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
