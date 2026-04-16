#ifndef PYTHONIC_INCLUDE_NUMPY_ISCOMPLEX_HPP
#define PYTHONIC_INCLUDE_NUMPY_ISCOMPLEX_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/traits.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace wrapper
  {
    template <class I>
    std::enable_if_t<types::is_complex<I>::value, bool> iscomplex(I const &a);

    template <class I>
    constexpr std::enable_if_t<!types::is_complex<I>::value, bool> iscomplex(I const &a);
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME iscomplex
#define NUMPY_NARY_FUNC_SYM wrapper::iscomplex
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
