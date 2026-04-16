#ifndef PYTHONIC_INCLUDE_NUMPY_HSTACK_HPP
#define PYTHONIC_INCLUDE_NUMPY_HSTACK_HPP

#include <pythonic/include/numpy/concatenate.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class ArraySequence>
  auto hstack(ArraySequence &&seq) -> decltype(concatenate(std::forward<ArraySequence>(seq), 1));

  DEFINE_FUNCTOR(pythonic::numpy, hstack);
} // namespace numpy
PYTHONIC_NS_END

#endif
