#ifndef PYTHONIC_INCLUDE_NUMPY_ISINF_HPP
#define PYTHONIC_INCLUDE_NUMPY_ISINF_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    bool isinf(T const &v);

    template <class T>
    bool isinf(std::complex<T> const &v);
  } // namespace wrapper
#define NUMPY_NARY_FUNC_NAME isinf
#define NUMPY_NARY_FUNC_SYM wrapper::isinf
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
