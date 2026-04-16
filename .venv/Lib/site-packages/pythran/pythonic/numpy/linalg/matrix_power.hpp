#ifndef PYTHONIC_NUMPY_LINALG_MATRIX_POWER_HPP
#define PYTHONIC_NUMPY_LINALG_MATRIX_POWER_HPP

#include "pythonic/include/numpy/linalg/matrix_power.hpp"

#include "pythonic/numpy/array.hpp"
#include "pythonic/numpy/asarray.hpp"
#include "pythonic/numpy/dot.hpp"
#include "pythonic/numpy/identity.hpp"

#include "pythonic/builtins/NotImplementedError.hpp"
PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace linalg
  {

    namespace details
    {

      template <class E>
      E fast_pow(E const &base, long n)
      {
        if (n == 1)
          return base;
        if (n == 2)
          return numpy::functor::dot{}(base, base);
        if (n == 3) {
          auto tmp = numpy::functor::dot{}(base, base);
          return numpy::functor::dot{}(tmp, base);
        }
        // starting from here, we know for sure that tmp will point to newly
        // allocated memory
        // this is used to optimize in-place dot computation in the odd case
        auto tmp = fast_pow(base, n / 2);
        if (n & 1) {
          auto next = numpy::functor::dot{}(tmp, tmp);
          return numpy::functor::dot{}(base, next, tmp);
        } else {
          return numpy::functor::dot{}(tmp, tmp);
        }
      }
    } // namespace details

    template <class E>
    auto matrix_power(E const &expr, long n) -> decltype(numpy::functor::array{}(expr))
    {
      if (n == 0)
        return numpy::functor::identity{}(expr.template shape<0>(),
                                          types::dtype_t<typename E::dtype>{});
      if (n > 0) {
        auto base = numpy::functor::asarray{}(expr);
        return details::fast_pow(base, n);
      }
      throw pythonic::builtins::NotImplementedError("negative power");
    }
  } // namespace linalg
} // namespace numpy
PYTHONIC_NS_END

#endif
