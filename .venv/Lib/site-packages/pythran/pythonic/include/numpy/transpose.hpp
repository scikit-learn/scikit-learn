#ifndef PYTHONIC_INCLUDE_NUMPY_TRANSPOSE_HPP
#define PYTHONIC_INCLUDE_NUMPY_TRANSPOSE_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/numpy_expr.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/nested_container.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E>
  types::numpy_texpr<types::broadcasted<E>> transpose(types::broadcasted<E> const &arr)
  {
    return {arr};
  }

  template <class E>
  E transpose(types::numpy_texpr<E> const &arr)
  {
    return arr.arg;
  }

  template <class E>
  std::enable_if_t<E::value == 2, types::numpy_texpr<E>> transpose(E const &arr)
  {
    return {arr};
  }

  template <class E>
  std::enable_if_t<E::value == 1, E> transpose(E const &arr)
  {
    return arr;
  }

  template <class T, class pS>
  std::enable_if_t<(std::tuple_size<pS>::value > 2),
                   types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>>
  transpose(types::ndarray<T, pS> const &a);

  template <class T, class pS, size_t M>
  types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>
  transpose(types::ndarray<T, pS> const &a, types::array_tuple<long, M> const &t);

  template <class T, class pS, class... Args>
  types::ndarray<T, types::array_tuple<long, 1 + sizeof...(Args)>>
  transpose(types::ndarray<T, pS> const &a, long index, Args const &...indices)
  {
    return transpose(a, types::array_tuple<long, 1 + sizeof...(Args)>{{index, (long)indices...}});
  }

  template <class T>
  struct _transpose;

  template <class Op, class... Args>
  auto transpose(types::numpy_expr<Op, Args...> const &expr)
      -> decltype(_transpose<types::numpy_expr<Op, Args...>>{}(
          expr, std::make_index_sequence<sizeof...(Args)>()))
  {
    return _transpose<types::numpy_expr<Op, Args...>>{}(
        expr, std::make_index_sequence<sizeof...(Args)>());
  }
  template <class Op, class... Args>
  struct _transpose<types::numpy_expr<Op, Args...>> {
    template <size_t... Is>
    auto operator()(types::numpy_expr<Op, Args...> const &expr, std::index_sequence<Is...>)
        -> decltype(Op{}(transpose(std::get<Is>(expr.args))...))
    {
      return Op{}(transpose(std::get<Is>(expr.args))...);
    }
  };

  template <class E>
  auto transpose(E const &expr) -> std::enable_if_t<
      (E::value > 2),
      decltype(transpose(types::ndarray<typename E::dtype, typename E::shape_t>{expr}))>
  {
    return transpose(types::ndarray<typename E::dtype, typename E::shape_t>{expr});
  }

  DEFINE_FUNCTOR(pythonic::numpy, transpose);
} // namespace numpy
PYTHONIC_NS_END

#endif
