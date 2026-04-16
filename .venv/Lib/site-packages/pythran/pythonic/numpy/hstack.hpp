#ifndef PYTHONIC_NUMPY_HSTACK_HPP
#define PYTHONIC_NUMPY_HSTACK_HPP

#include <pythonic/include/numpy/hstack.hpp>
#include <pythonic/numpy/concatenate.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class ArraySequence>
  auto hstack(ArraySequence &&seq) -> decltype(concatenate(std::forward<ArraySequence>(seq), 1))
  {
    auto constexpr concatenate_axis =
        (decltype(concatenate(std::forward<ArraySequence>(seq), 1))::value != 1);
    return concatenate(std::forward<ArraySequence>(seq), concatenate_axis);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
