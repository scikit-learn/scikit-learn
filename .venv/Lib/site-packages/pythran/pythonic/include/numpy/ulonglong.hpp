#ifndef PYTHONIC_INCLUDE_NUMPY_ULONGLONG_HPP
#define PYTHONIC_INCLUDE_NUMPY_ULONGLONG_HPP

#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    unsigned long long ulonglong();
    template <class V>
    unsigned long long ulonglong(V v);
  } // namespace details

#define NUMPY_NARY_FUNC_NAME ulonglong
#define NUMPY_NARY_FUNC_SYM details::ulonglong
#define NUMPY_NARY_EXTRA_METHOD using type = unsigned long long;
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
