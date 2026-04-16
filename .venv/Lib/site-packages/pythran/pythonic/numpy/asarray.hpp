#ifndef PYTHONIC_NUMPY_ASARRAY_HPP
#define PYTHONIC_NUMPY_ASARRAY_HPP

#include "pythonic/include/numpy/asarray.hpp"

#include "pythonic/numpy/array.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E, class dtype>
  template <class... Types>
  auto _asarray<E, dtype>::operator()(Types &&...args)
      -> decltype(array(std::forward<Types>(args)...))
  {
    return array(std::forward<Types>(args)...);
  }

  template <class T, class pS>
  template <class F, class dtype>
  F &&_asarray<types::ndarray<T, pS>, T>::operator()(F &&a, dtype)
  {
    return std::forward<F>(a);
  }

  template <class E>
  auto asarray(E &&e, types::none_type d)
      -> decltype(_asarray<std::decay_t<E>, typename types::dtype_of<std::decay_t<E>>::type>{}(
          std::forward<E>(e)))
  {
    return _asarray<std::decay_t<E>, typename types::dtype_of<std::decay_t<E>>::type>{}(
        std::forward<E>(e));
  }

  template <class E, class dtype>
  auto asarray(E &&e, dtype d)
      -> decltype(_asarray<std::decay_t<E>, typename dtype::type>{}(std::forward<E>(e), d))
  {
    return _asarray<std::decay_t<E>, typename dtype::type>{}(std::forward<E>(e), d);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
