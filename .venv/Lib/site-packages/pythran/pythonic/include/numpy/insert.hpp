#ifndef PYTHONIC_INCLUDE_NUMPY_INSERT_HPP
#define PYTHONIC_INCLUDE_NUMPY_INSERT_HPP

#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/traits.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class T, class pS, class I, class F>
  std::enable_if_t<types::is_iterable<I>::value && types::is_iterable<F>::value,
                   types::ndarray<T, types::pshape<long>>>
  insert(types::ndarray<T, pS> in, I const &indices, F const &data,
         types::none_type axis = builtins::None);

  template <class T, class pS, class I, class F>
  std::enable_if_t<types::is_iterable<I>::value && !types::is_iterable<F>::value,
                   types::ndarray<T, types::pshape<long>>>
  insert(types::ndarray<T, pS> in, I const &indices, F const &data,
         types::none_type axis = builtins::None);

  template <class T, class pS, class I, class F>
  std::enable_if_t<!types::is_iterable<I>::value && types::is_iterable<F>::value,
                   types::ndarray<T, types::pshape<long>>>
  insert(types::ndarray<T, pS> in, I const &indices, F const &data,
         types::none_type axis = builtins::None);

  template <class T, class pS, class I, class F>
  std::enable_if_t<!types::is_iterable<I>::value && !types::is_iterable<F>::value,
                   types::ndarray<T, types::pshape<long>>>
  insert(types::ndarray<T, pS> in, I const &indices, F const &data,
         types::none_type axis = builtins::None);

  template <class E, class... Args>
  E insert(E, Args const &...);

  DEFINE_FUNCTOR(pythonic::numpy, insert);
} // namespace numpy
PYTHONIC_NS_END

#endif
