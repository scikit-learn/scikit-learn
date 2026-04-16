#ifndef PYTHONIC_INCLUDE_NUMPY_FILL_DIAGONAL_HPP
#define PYTHONIC_INCLUDE_NUMPY_FILL_DIAGONAL_HPP

#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::none_type fill_diagonal(E &&, typename std::decay_t<E>::dtype);

  DEFINE_FUNCTOR(pythonic::numpy, fill_diagonal)
} // namespace numpy
PYTHONIC_NS_END

#endif
