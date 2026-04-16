#ifndef PYTHONIC_INCLUDE_NUMPY_ASSCALAR_HPP
#define PYTHONIC_INCLUDE_NUMPY_ASSCALAR_HPP

#include "pythonic/include/numpy/asarray.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  using asscalar_result_type = std::conditional_t<
      std::is_integral<T>::value, long,
      std::conditional_t<std::is_floating_point<T>::value, double, std::complex<double>>>;

  template <class E>
  asscalar_result_type<typename E::dtype> asscalar(E const &expr);

  DEFINE_FUNCTOR(pythonic::numpy, asscalar);
} // namespace numpy
PYTHONIC_NS_END

#endif
