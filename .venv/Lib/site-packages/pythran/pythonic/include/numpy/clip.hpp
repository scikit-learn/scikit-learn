#ifndef PYTHONIC_INCLUDE_NUMPY_CLIP_HPP
#define PYTHONIC_INCLUDE_NUMPY_CLIP_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <xsimd/xsimd.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{

  // special private handler for clip(v, None, m) {

#define NUMPY_NARY_FUNC_NAME _clip_max
#define NUMPY_NARY_FUNC_SYM xsimd::min
#include "pythonic/include/types/numpy_nary_expr.hpp"

  // }

  namespace wrapper
  {
    template <class T, class Mi, class Ma>
    typename __combined<T, Mi, Ma>::type clip(T const &v, Mi a_min, Ma a_max);

    template <class T, class Mi>
    typename __combined<T, Mi>::type clip(T const &v, Mi a_min);
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME clip
#define NUMPY_NARY_FUNC_SYM wrapper::clip
#define NUMPY_NARY_EXTRA_METHOD                                                                    \
  template <typename T, class Mi>                                                                  \
  auto operator()(T &&v, Mi &&a_min, types::none_type)                                             \
      ->decltype((*this)(std::forward<T>(v), std::forward<Mi>(a_min)))                             \
  {                                                                                                \
    return (*this)(std::forward<T>(v), std::forward<Mi>(a_min));                                   \
  }                                                                                                \
  template <typename T, class Ma>                                                                  \
  auto operator()(T &&v, types::none_type, Ma &&a_max)                                             \
      ->decltype(_clip_max{}(std::forward<T>(v), std::forward<Ma>(a_max)))                         \
  {                                                                                                \
    return _clip_max{}(std::forward<T>(v), std::forward<Ma>(a_max));                               \
  }
#include "pythonic/include/types/numpy_nary_expr.hpp"

} // namespace numpy
PYTHONIC_NS_END

#endif
