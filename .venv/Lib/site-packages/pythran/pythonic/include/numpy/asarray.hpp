#ifndef PYTHONIC_INCLUDE_NUMPY_ASARRAY_HPP
#define PYTHONIC_INCLUDE_NUMPY_ASARRAY_HPP

#include "pythonic/include/numpy/array.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class dtype>
  struct _asarray {
    template <class... Types>
    auto operator()(Types &&...args) -> decltype(array(std::forward<Types>(args)...));
  };

  template <class T, class pS>
  struct _asarray<types::ndarray<T, pS>, T> {
    template <class F, class dtype = types::none_type>
    F &&operator()(F &&a, dtype d = types::none_type());
  };

  template <class E>
  auto asarray(E &&e, types::none_type d = types::none_type())
      -> decltype(_asarray<std::decay_t<E>, typename types::dtype_of<std::decay_t<E>>::type>{}(
          std::forward<E>(e)));

  template <class E, class dtype>
  auto asarray(E &&e, dtype d)
      -> decltype(_asarray<std::decay_t<E>, typename dtype::type>{}(std::forward<E>(e), d));

  DEFINE_FUNCTOR(pythonic::numpy, asarray);
} // namespace numpy
PYTHONIC_NS_END

#endif
