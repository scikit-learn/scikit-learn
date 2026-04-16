#ifndef PYTHONIC_NUMPY_NDINDEX_HPP
#define PYTHONIC_NUMPY_NDINDEX_HPP

#include "pythonic/include/numpy/ndindex.hpp"

#include "pythonic/utils/functor.hpp"

#include <numeric>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <size_t N>
  ndindex_iterator<N>::ndindex_iterator()
  {
  }

  template <size_t N>
  ndindex_iterator<N>::ndindex_iterator(types::array_tuple<long, N> const &shape, long first)
      : index(first), shape(shape)
  {
  }

  template <size_t N>
  types::array_tuple<long, N> ndindex_iterator<N>::operator*() const
  {
    types::array_tuple<long, N> out;
    long mult = 1;
    for (long j = N - 1; j > 0; j--) {
      out[j] = (index / mult) % shape[j];
      mult *= shape[j];
    }
    out[0] = index / mult;
    return out;
  }

  template <size_t N>
  ndindex_iterator<N> &ndindex_iterator<N>::operator++()
  {
    ++index;
    return *this;
  }

  template <size_t N>
  ndindex_iterator<N> &ndindex_iterator<N>::operator+=(long n)
  {
    index += n;
    return *this;
  }

  template <size_t N>
  bool ndindex_iterator<N>::operator!=(ndindex_iterator<N> const &other) const
  {
    return index != other.index;
  }

  template <size_t N>
  bool ndindex_iterator<N>::operator<(ndindex_iterator<N> const &other) const
  {
    return index < other.index;
  }

  template <size_t N>
  long ndindex_iterator<N>::operator-(ndindex_iterator<N> const &other) const
  {
    return index - other.index;
  }

  template <size_t N>
  _ndindex<N>::_ndindex()
  {
  }

  template <size_t N>
  _ndindex<N>::_ndindex(types::array_tuple<long, N> const &shape)
      : ndindex_iterator<N>(shape, 0), shape(shape),
        end_iter(shape, std::accumulate(shape.begin(), shape.end(), 1L, std::multiplies<long>()))
  {
  }

  template <size_t N>
  typename _ndindex<N>::iterator &_ndindex<N>::begin()
  {
    return *this;
  }

  template <size_t N>
  typename _ndindex<N>::iterator const &_ndindex<N>::begin() const
  {
    return *this;
  }

  template <size_t N>
  typename _ndindex<N>::iterator _ndindex<N>::end() const
  {
    return end_iter;
  }

  template <class... Types>
  _ndindex<sizeof...(Types)> ndindex(Types... args)
  {
    return {types::make_tuple(args...)};
  }

  template <size_t N>
  _ndindex<N> ndindex(types::array_tuple<long, N> const &args)
  {
    return {args};
  }
  template <class... Tys>
  _ndindex<sizeof...(Tys)> ndindex(types::pshape<Tys...> const &args)
  {
    return {args};
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
