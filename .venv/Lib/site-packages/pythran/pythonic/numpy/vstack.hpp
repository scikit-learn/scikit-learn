#ifndef PYTHONIC_NUMPY_VSTACK_HPP
#define PYTHONIC_NUMPY_VSTACK_HPP

#include <pythonic/include/numpy/vstack.hpp>
#include <pythonic/numpy/concatenate.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class ArraySequence>
  auto vstack(ArraySequence &&seq)
      -> std::enable_if_t<(impl::vstack_helper<ArraySequence>::value > 1),
                          impl::vstack_helper<ArraySequence>>
  {

    return concatenate(std::forward<ArraySequence>(seq), 0);
  }

  template <class ArraySequence>
  auto vstack(ArraySequence &&seq)
      -> std::enable_if_t<(impl::vstack_helper<ArraySequence>::value == 1),
                          decltype(std::declval<impl::vstack_helper<ArraySequence>>().reshape(
                              std::declval<types::array_tuple<long, 2>>()))>
  {
    auto &&temp = concatenate(std::forward<ArraySequence>(seq), 0);
    long const seq_size = seq.size(), temp_size = temp.size();
    types::array_tuple<long, 2> new_shape{{seq_size, temp_size / seq_size}};
    return temp.reshape(new_shape);
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
