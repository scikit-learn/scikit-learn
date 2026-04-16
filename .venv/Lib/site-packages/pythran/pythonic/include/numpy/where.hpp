#ifndef PYTHONIC_INCLUDE_NUMPY_WHERE_HPP
#define PYTHONIC_INCLUDE_NUMPY_WHERE_HPP

#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/numpy/copy.hpp"
#include "pythonic/include/numpy/nonzero.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace impl
  {
    template <class E, class F, class G>
    typename __combined<F, G>::type where(E const &cond, F const &true_, G const &false_);
  }

#define NUMPY_NARY_EXTRA_METHOD                                                                    \
  template <class E>                                                                               \
  auto operator()(E &&expr)->decltype(nonzero{}(std::forward<E>(expr)))                            \
  {                                                                                                \
    return nonzero{}(std::forward<E>(expr));                                                       \
  }

#define NUMPY_NARY_FUNC_NAME where
#define NUMPY_NARY_FUNC_SYM impl::where
#define NUMPY_NARY_RESHAPE_MODE reshape_type
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
