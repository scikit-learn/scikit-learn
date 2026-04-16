#ifndef PYTHONIC_NUMPY_INTP_HPP
#define PYTHONIC_NUMPY_INTP_HPP

#include "pythonic/include/numpy/intp.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline intptr_t intp()
    {
      return intptr_t();
    }

    template <class V>
    intptr_t intp(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME intp
#define NUMPY_NARY_FUNC_SYM details::intp
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
