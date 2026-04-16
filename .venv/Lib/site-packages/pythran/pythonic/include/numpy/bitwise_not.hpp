#ifndef PYTHONIC_INCLUDE_NUMPY_BITWISENOT_HPP
#define PYTHONIC_INCLUDE_NUMPY_BITWISENOT_HPP

#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{

  namespace wrapper
  {
    template <class A>
    decltype(~std::declval<A const &>()) bitwise_not(A const &a);

    bool bitwise_not(bool t0);
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME bitwise_not
#define NUMPY_NARY_FUNC_SYM wrapper::bitwise_not
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

// ndarray have to be include after as bitwise_not is used as a numpy_operator
#include "pythonic/include/types/ndarray.hpp"

#endif
