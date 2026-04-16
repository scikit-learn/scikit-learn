#ifndef PYTHONIC_NUMPY_ISNAN_HPP
#define PYTHONIC_NUMPY_ISNAN_HPP

#include "pythonic/include/numpy/isnan.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    bool isnan(std::complex<T> const &v)
    {
      return std::isnan(v.real()) || std::isnan(v.imag());
    }

    template <class T>
    auto isnan(T const &v) -> std::enable_if_t<std::is_floating_point<std::decay_t<T>>::value, bool>
    {
      return std::isnan(v);
    }

    template <class T>
    auto isnan(T const &v)
        -> std::enable_if_t<!std::is_floating_point<std::decay_t<T>>::value, bool>
    {
      return false;
    }
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME isnan
#define NUMPY_NARY_FUNC_SYM wrapper::isnan
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
