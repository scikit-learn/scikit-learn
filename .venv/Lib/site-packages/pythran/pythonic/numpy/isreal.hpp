#ifndef PYTHONIC_NUMPY_ISREAL_HPP
#define PYTHONIC_NUMPY_ISREAL_HPP

#include "pythonic/include/numpy/isreal.hpp"

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
    std::enable_if_t<types::is_complex<I>::value, bool> isreal(I const &a)
    {
      return a.imag() == 0.;
    }

    template <class I>
    std::enable_if_t<!types::is_complex<I>::value, bool> isreal(I const &a)
    {
      return true;
    }
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME isreal
#define NUMPY_NARY_FUNC_SYM wrapper::isreal
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
