#ifndef PYTHONIC_INCLUDE_NUMPY_FIX_HPP
#define PYTHONIC_INCLUDE_NUMPY_FIX_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class E>
    E fix(E const &e)
    {
      if (std::is_integral<E>::value)
        return e;
      else
        return std::trunc(e);
    }
  } // namespace wrapper
#define NUMPY_NARY_FUNC_NAME fix
#define NUMPY_NARY_FUNC_SYM wrapper::fix
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
