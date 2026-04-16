#ifndef PYTHONIC_NUMPY_ANY_HPP
#define PYTHONIC_NUMPY_ANY_HPP

#include "pythonic/include/numpy/any.hpp"

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/numpy/add.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  bool _any(E const &e, utils::int_<1>)
  {
    return std::any_of(e.begin(), e.end(), [](typename E::dtype elt) -> bool { return elt; });
  }

  template <class E, size_t N>
  bool _any(E const &e, utils::int_<N>)
  {
    for (auto &&elt : e)
      if (_any(elt, utils::int_<N - 1>())) {
        return true;
      }
    return false;
  }

  template <class E>
  std::enable_if_t<types::is_numexpr_arg<E>::value, bool> any(E const &expr, types::none_type)
  {
    return _any(expr, utils::int_<E::value>());
  }

  template <class E>
  std::enable_if_t<std::is_scalar<E>::value || types::is_complex<E>::value, bool>
  any(E const &expr, types::none_type)
  {
    return expr;
  }

  template <class E>
  auto any(E const &array, long axis)
      -> std::enable_if_t<std::is_scalar<E>::value || types::is_complex<E>::value,
                          decltype(any(array))>
  {
    if (axis != 0)
      throw types::ValueError("axis out of bounds");
    return any(array);
  }

  template <class E>
  auto any(E const &array, long axis) -> std::enable_if_t<E::value == 1, decltype(any(array))>
  {
    if (axis != 0)
      throw types::ValueError("axis out of bounds");
    return any(array);
  }

  template <class E>
  std::enable_if_t<E::value != 1,
                   types::ndarray<typename E::dtype, types::array_tuple<long, E::value - 1>>>
  any(E const &array, long axis)
  {
    constexpr long N = E::value;
    using T = typename E::dtype;
    if (axis < 0 || axis >= long(N))
      throw types::ValueError("axis out of bounds");
    if (axis == 0) {
      types::array_tuple<long, N> shp;
      shp[0] = 1;
      sutils::copy_shape<1, 0>(shp, array, std::make_index_sequence<N - 1>());
      types::ndarray<bool, types::array_tuple<long, N>> out(shp, false);
      return std::accumulate(array.begin(), array.end(), *out.begin(), numpy::functor::add());
    } else {
      types::array_tuple<long, N - 1> shp;
      sutils::copy_shape<0, 0>(shp, array, std::make_index_sequence<N - 1>());
      types::ndarray<bool, types::array_tuple<long, N - 1>> anyy(shp, builtins::None);
      std::transform(array.begin(), array.end(), anyy.begin(),
                     [=](types::ndarray<T, types::array_tuple<long, N - 1>> const &other) {
                       return any(other, axis - 1);
                     });
      return anyy;
    }
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
