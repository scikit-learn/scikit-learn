#ifndef PYTHONIC_INCLUDE_NUMPY_UNRAVEL_INDEX_HPP
#define PYTHONIC_INCLUDE_NUMPY_UNRAVEL_INDEX_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/str.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class S>
  std::enable_if_t<std::is_scalar<E>::value, types::array_tuple<long, std::tuple_size<S>::value>>
  unravel_index(E const &expr, S const &shape, types::str const &order = "C");

  DEFINE_FUNCTOR(pythonic::numpy, unravel_index);
} // namespace numpy
PYTHONIC_NS_END

#endif
