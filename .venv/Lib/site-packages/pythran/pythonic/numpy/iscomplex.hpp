#ifndef PYTHONIC_NUMPY_ISCOMPLEX_HPP
#define PYTHONIC_NUMPY_ISCOMPLEX_HPP

#include "pythonic/include/numpy/iscomplex.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/traits.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace wrapper
  {
    template <class I>
    std::enable_if_t<types::is_complex<I>::value, bool> iscomplex(I const &a)
    {
      return a.imag() != 0.;
    }

    template <class I>
    constexpr std::enable_if_t<!types::is_complex<I>::value, bool> iscomplex(I const &a)
    {
      return false;
    }
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME iscomplex
#define NUMPY_NARY_FUNC_SYM wrapper::iscomplex
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
