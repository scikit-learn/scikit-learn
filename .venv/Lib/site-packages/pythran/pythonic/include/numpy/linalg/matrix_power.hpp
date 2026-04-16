#ifndef PYTHONIC_INCLUDE_NUMPY_LINALG_MATRIX_POWER_HPP
#define PYTHONIC_INCLUDE_NUMPY_LINALG_MATRIX_POWER_HPP

#include "pythonic/include/numpy/array.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace linalg
  {
    template <class E>
    auto matrix_power(E const &expr, long n) -> decltype(numpy::functor::array{}(expr));

    DEFINE_FUNCTOR(pythonic::numpy::linalg, matrix_power);
  } // namespace linalg
} // namespace numpy
PYTHONIC_NS_END

#endif
