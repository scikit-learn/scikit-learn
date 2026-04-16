#ifndef PYTHONIC_INCLUDE_NUMPY_ARRAYSPLIT_HPP
#define PYTHONIC_INCLUDE_NUMPY_ARRAYSPLIT_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::list<
      typename assignable<decltype(std::declval<E>()[types::fast_contiguous_slice()])>::type>
  array_split(E const &a, long nb_split);

  template <class E, class I>
  std::enable_if_t<types::is_iterable<I>::value,
                   types::list<typename assignable<
                       decltype(std::declval<E>()[types::fast_contiguous_slice()])>::type>>
  array_split(E const &a, I const &split_mask);

  NUMPY_EXPR_TO_NDARRAY0_DECL(array_split);
  DEFINE_FUNCTOR(pythonic::numpy, array_split);
} // namespace numpy
PYTHONIC_NS_END

#endif
