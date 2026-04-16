#ifndef PYTHONIC_INCLUDE_NUMPY_SQUARE_HPP
#define PYTHONIC_INCLUDE_NUMPY_SQUARE_HPP

#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

#include <complex>

namespace wrapper
{
}

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace wrapper
  {
    template <class T>
    auto square(T const &arg) -> decltype(arg * arg)
    {
      return arg * arg;
    }
    template <class T>
    std::complex<T> square(std::complex<T> const &arg)
    {
      T r = arg.real(), i = arg.imag();
      auto t = r * i;
      auto r2 = r * r;
      auto i2 = i * i;
      return {r2 - i2, t + t};
    }
  } // namespace wrapper

#define NUMPY_NARY_FUNC_NAME square
#define NUMPY_NARY_FUNC_SYM wrapper::square
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#endif
