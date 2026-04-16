#ifndef PYTHONIC_NUMPY_ULONGLONG_HPP
#define PYTHONIC_NUMPY_ULONGLONG_HPP

#include "pythonic/include/numpy/ulonglong.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline unsigned long long ulonglong()
    {
      return {};
    }

    template <class V>
    unsigned long long ulonglong(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME ulonglong
#define NUMPY_NARY_FUNC_SYM details::ulonglong
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
