#ifndef PYTHONIC_INCLUDE_NUMPY_VSTACK_HPP
#define PYTHONIC_INCLUDE_NUMPY_VSTACK_HPP

#include <pythonic/include/numpy/concatenate.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace impl
  {
    template <class T>
    using vstack_helper = decltype(concatenate(std::declval<T>(), 0));
  }

  template <class ArraySequence>
  auto vstack(ArraySequence &&seq)
      -> std::enable_if_t<(impl::vstack_helper<ArraySequence>::value > 1),
                          impl::vstack_helper<ArraySequence>>;

  // according to the numpy.vstack doc:
  //  Equivalent to ``np.concatenate(tup, axis=0)`` if `tup` contains arrays
  //  that
  //      are at least 2-dimensional.
  //
  // the enable if is there to match this behavior
  template <class ArraySequence>
  auto vstack(ArraySequence &&seq)
      -> std::enable_if_t<(impl::vstack_helper<ArraySequence>::value == 1),
                          decltype(std::declval<impl::vstack_helper<ArraySequence>>().reshape(
                              std::declval<types::array_tuple<long, 2>>()))>;

  DEFINE_FUNCTOR(pythonic::numpy, vstack);
} // namespace numpy
PYTHONIC_NS_END

#endif
