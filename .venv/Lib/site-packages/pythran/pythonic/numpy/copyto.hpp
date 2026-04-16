#ifndef PYTHONIC_NUMPY_COPYTO_HPP
#define PYTHONIC_NUMPY_COPYTO_HPP

#include "pythonic//numpy/asarray.hpp"
#include "pythonic/include/numpy/copyto.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  template <class T, class pS, class E>
  types::none_type copyto(types::ndarray<T, pS> &out, E const &expr)
  {
    using out_type = types::ndarray<T, pS>;
    if (may_overlap(out, expr)) {
      auto aexpr = asarray(expr);
      utils::broadcast_copy<
          out_type &, decltype(aexpr), out_type::value,
          (int)out_type::value - (int)utils::dim_of<E>::value,
          out_type::is_vectorizable &&
              std::is_same<typename out_type::dtype, typename types::dtype_of<E>::type>::value &&
              types::is_vectorizable<E>::value>(out, aexpr);
    } else {
      utils::broadcast_copy<
          out_type &, E, out_type::value, (int)out_type::value - (int)utils::dim_of<E>::value,
          out_type::is_vectorizable &&
              std::is_same<typename out_type::dtype, typename types::dtype_of<E>::type>::value &&
              types::is_vectorizable<E>::value>(out, expr);
    }
    return {};
  }

  template <class T, class pS, class E>
  types::none_type copyto(types::ndarray<T, pS> &&out, E const &expr)
  {
    return copyto(out, expr);
  }

  template <class T, class pS, class E>
  types::none_type copyto(types::numpy_texpr<types::ndarray<T, pS>> &out, E const &expr)
  {
    using out_type = types::numpy_texpr<types::ndarray<T, pS>>;
    if (may_overlap(out, expr)) {
      auto aexpr = asarray(expr);
      utils::broadcast_copy<
          out_type &, decltype(aexpr), out_type::value,
          (int)out_type::value - (int)utils::dim_of<E>::value,
          out_type::is_vectorizable &&
              std::is_same<typename out_type::dtype, typename types::dtype_of<E>::type>::value &&
              types::is_vectorizable<E>::value>(out, aexpr);
    } else {
      utils::broadcast_copy<
          out_type &, E, out_type::value, (int)out_type::value - (int)utils::dim_of<E>::value,
          out_type::is_vectorizable &&
              std::is_same<typename out_type::dtype, typename types::dtype_of<E>::type>::value &&
              types::is_vectorizable<E>::value>(out, expr);
    }
    return {};
  }

  template <class T, class pS, class E>
  types::none_type copyto(types::numpy_texpr<types::ndarray<T, pS>> &&out, E const &expr)
  {
    return copyto(out, expr);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
