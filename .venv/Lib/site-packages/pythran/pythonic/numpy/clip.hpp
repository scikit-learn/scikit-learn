#ifndef PYTHONIC_NUMPY_CLIP_HPP
#define PYTHONIC_NUMPY_CLIP_HPP

#include "pythonic/include/numpy/clip.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME _clip_max
#define NUMPY_NARY_FUNC_SYM xsimd::min
#include "pythonic/types/numpy_nary_expr.hpp"

  namespace wrapper
  {
    template <class T, class Mi, class Ma>
    typename __combined<T, Mi, Ma>::type clip(T const &v, Mi a_min, Ma a_max)
    {
      if (v < a_min)
        return a_min;
      else if (v > a_max)
        return a_max;
      else
        return v;
    }

    template <class T, class Mi>
    typename __combined<T, Mi>::type clip(T const &v, Mi a_min)
    {
      if (v < a_min)
        return a_min;
      else
        return v;
    }
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME clip
#define NUMPY_NARY_FUNC_SYM wrapper::clip
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
