#ifndef PYTHONIC_INCLUDE_NUMPY_SPLIT_HPP
#define PYTHONIC_INCLUDE_NUMPY_SPLIT_HPP

#include "pythonic/include/numpy/array_split.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::list<types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>>
  split(types::ndarray<T, pS> const &a, long nb_split);

  template <class T, class pS, class I>
  std::enable_if_t<
      types::is_iterable<I>::value,
      types::list<types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>>>
  split(types::ndarray<T, pS> const &a, I const &split_mask);

  template <class E, class I>
  types::list<types::ndarray<typename E::dtype, types::array_tuple<long, E::value>>>
  split(E const &a, I const &);

  DEFINE_FUNCTOR(pythonic::numpy, split);
} // namespace numpy
PYTHONIC_NS_END

#endif
