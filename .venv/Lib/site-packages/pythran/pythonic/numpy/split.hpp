#ifndef PYTHONIC_NUMPY_SPLIT_HPP
#define PYTHONIC_NUMPY_SPLIT_HPP

#include "pythonic/include/numpy/split.hpp"

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/numpy/array_split.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::list<types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>>
  split(types::ndarray<T, pS> const &a, long nb_split)
  {
    if (a.flat_size() % nb_split != 0)
      throw types::ValueError("array split does ! result in an equal division");
    return array_split(a, nb_split);
  }

  template <class T, class pS, class I>
  std::enable_if_t<
      types::is_iterable<I>::value,
      types::list<types::ndarray<T, types::array_tuple<long, std::tuple_size<pS>::value>>>>
  split(types::ndarray<T, pS> const &a, I const &split_mask)
  {
    return array_split(a, split_mask);
  }

  template <class E, class I>
  types::list<types::ndarray<typename E::dtype, types::array_tuple<long, E::value>>>
  split(E const &a, I const &)
  {
    throw std::runtime_error("split only partially implemented");
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
