#ifndef PYTHONIC_NUMPY_BOOL_HPP
#define PYTHONIC_NUMPY_BOOL_HPP

#include "pythonic/include/numpy/bool_.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline bool bool_()
    {
      return bool();
    }

    template <class V>
    bool bool_(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME bool_
#define NUMPY_NARY_FUNC_SYM details::bool_
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
