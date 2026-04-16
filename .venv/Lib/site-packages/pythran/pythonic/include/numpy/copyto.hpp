#ifndef PYTHONIC_INCLUDE_NUMPY_COPYTO_HPP
#define PYTHONIC_INCLUDE_NUMPY_COPYTO_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  template <class T, class pS, class E>
  types::none_type copyto(types::ndarray<T, pS> &out, E const &expr);

  template <class T, class pS, class E>
  types::none_type copyto(types::ndarray<T, pS> &&out, E const &expr);

  template <class T, class pS, class E>
  types::none_type copyto(types::numpy_texpr<types::ndarray<T, pS>> &out, E const &expr);

  template <class T, class pS, class E>
  types::none_type copyto(types::numpy_texpr<types::ndarray<T, pS>> &&out, E const &expr);

  // pythran extensions
  template <class E, class F>
  types::none_type copyto(E &out, F const &expr)
  {
    out[types::fast_contiguous_slice(0, types::none_type{})] = expr;
    return {};
  }

  DEFINE_FUNCTOR(pythonic::numpy, copyto);
} // namespace numpy
PYTHONIC_NS_END

#endif
