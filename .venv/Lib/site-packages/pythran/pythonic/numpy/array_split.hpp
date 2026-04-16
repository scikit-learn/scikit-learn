#ifndef PYTHONIC_NUMPY_ARRAYSPLIT_HPP
#define PYTHONIC_NUMPY_ARRAYSPLIT_HPP

#include "pythonic/include/numpy/array_split.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::list<
      typename assignable<decltype(std::declval<E>()[types::fast_contiguous_slice()])>::type>
  array_split(E const &a, long nb_split)
  {
    long sz = a.template shape<0>();
    long n = (sz + nb_split - 1) / nb_split;
    long end = n * nb_split;
    long nb_full_split = nb_split;
    if (end != sz)
      nb_full_split -= (end - sz);
    types::list<
        typename assignable<decltype(std::declval<E>()[types::fast_contiguous_slice()])>::type>
        out(nb_split);

    long index = 0;
    for (long i = 0; i < nb_full_split; ++i, index += n)
      out[i] = a[types::fast_contiguous_slice(index, index + n)];
    for (long i = nb_full_split; i < nb_split; ++i, index += (n - 1))
      out[i] = a[types::fast_contiguous_slice(index, index + n - 1)];

    return out;
  }

  template <class E, class I>
  std::enable_if_t<types::is_iterable<I>::value,
                   types::list<typename assignable<
                       decltype(std::declval<E>()[types::fast_contiguous_slice()])>::type>>
  array_split(E const &a, I const &split_mask)
  {
    long sz = a.template shape<0>();
    types::list<
        typename assignable<decltype(std::declval<E>()[types::fast_contiguous_slice()])>::type>
        out(1 + split_mask.flat_size());
    long index = 0;
    auto inserter = out.begin();
    for (auto next_index : split_mask) {
      *inserter++ = a[types::fast_contiguous_slice(index, next_index)];
      index = next_index;
    }
    *inserter = a[types::fast_contiguous_slice(index, sz)];
    return out;
  }

  NUMPY_EXPR_TO_NDARRAY0_IMPL(array_split);
} // namespace numpy
PYTHONIC_NS_END

#endif
