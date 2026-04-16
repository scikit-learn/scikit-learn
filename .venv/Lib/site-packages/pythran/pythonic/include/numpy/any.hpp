#ifndef PYTHONIC_INCLUDE_NUMPY_ANY_HPP
#define PYTHONIC_INCLUDE_NUMPY_ANY_HPP

#include "pythonic/include/numpy/add.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  std::enable_if_t<types::is_numexpr_arg<E>::value, bool>
  any(E const &expr, types::none_type _ = types::none_type());

  template <class E>
  std::enable_if_t<std::is_scalar<E>::value || types::is_complex<E>::value, bool>
  any(E const &expr, types::none_type _ = types::none_type());

  template <class E>
  auto any(E const &array, long axis)
      -> std::enable_if_t<std::is_scalar<E>::value || types::is_complex<E>::value,
                          decltype(any(array))>;

  template <class E>
  auto any(E const &array, long axis) -> std::enable_if_t<E::value == 1, decltype(any(array))>;

  template <class E>
  std::enable_if_t<E::value != 1,
                   types::ndarray<typename E::dtype, types::array_tuple<long, E::value - 1>>>
  any(E const &array, long axis);

  DEFINE_FUNCTOR(pythonic::numpy, any);
} // namespace numpy
PYTHONIC_NS_END

#endif
