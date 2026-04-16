#ifndef PYTHONIC_INCLUDE_NUMPY_APPEND_HPP
#define PYTHONIC_INCLUDE_NUMPY_APPEND_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS, class F>
  std::enable_if_t<!types::is_dtype<F>::value,
                   types::ndarray<typename __combined<T, typename types::dtype_of<F>::type>::type,
                                  types::pshape<long>>>
  append(types::ndarray<T, pS> const &nto, F const &data);

  template <class T, class pS, class F>
  std::enable_if_t<types::is_dtype<F>::value,
                   types::ndarray<typename __combined<T, typename types::dtype_of<F>::type>::type,
                                  types::pshape<long>>>
  append(types::ndarray<T, pS> const &nto, F const &data);

  template <class T, class F>
  types::ndarray<typename __combined<typename types::dtype_of<T>::type,
                                     typename types::dtype_of<F>::type>::type,
                 types::pshape<long>>
  append(T const &to, F const &data);

  DEFINE_FUNCTOR(pythonic::numpy, append);
} // namespace numpy
PYTHONIC_NS_END

#endif
