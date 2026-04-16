#ifndef PYTHONIC_INCLUDE_NUMPY_PUT_HPP
#define PYTHONIC_INCLUDE_NUMPY_PUT_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class F, class T, class pS, class E>
  std::enable_if_t<types::is_numexpr_arg<F>::value, types::none_type>
  put(types::ndarray<T, pS> &expr, F const &ind, E const &v);

  template <class T, class pS>
  types::none_type put(types::ndarray<T, pS> &expr, long int ind, T const &v);

  template <class E, class M, class V>
  types::none_type put(E &, M const &, V const &);

  DEFINE_FUNCTOR(pythonic::numpy, put);
} // namespace numpy
PYTHONIC_NS_END

#endif
