#ifndef PYTHONIC_INCLUDE_NUMPY_SEARCHSORTED_HPP
#define PYTHONIC_INCLUDE_NUMPY_SEARCHSORTED_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/int_.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class T, class U>
  std::enable_if_t<!types::is_numexpr_arg<T>::value, long>
  searchsorted(U const &a, T const &v, types::str const &side = "left");

  template <class E, class T>
  std::enable_if_t<types::is_numexpr_arg<E>::value,
                   types::ndarray<long, types::array_tuple<long, E::value>>>
  searchsorted(T const &a, E const &v, types::str const &side = "left");

  DEFINE_FUNCTOR(pythonic::numpy, searchsorted);
} // namespace numpy
PYTHONIC_NS_END

#endif
