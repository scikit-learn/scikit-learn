#ifndef PYTHONIC_ITERTOOLS_PERMUTATIONS_HPP
#define PYTHONIC_ITERTOOLS_PERMUTATIONS_HPP

#include "pythonic/builtins/range.hpp"
#include "pythonic/include/itertools/permutations.hpp"
#include "pythonic/types/dynamic_tuple.hpp"
#include "pythonic/utils/functor.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace itertools
{

  template <class T, class H>
  permutations_iterator<T, H>::permutations_iterator()
  {
  }

  template <class T, class H>
  permutations_iterator<T, H>::permutations_iterator(pool_type const &iter, size_t num_elts,
                                                     bool end)
      : pool(iter), curr_permut(pool.size()), _size(num_elts), end(end)
  {
    std::iota(curr_permut.begin(), curr_permut.end(), 0);
    if (num_elts > iter.size()) {
      end = true;
    }
  }

  template <class T>
  types::dynamic_tuple<T> init_permut_from(size_t n, types::dynamic_tuple<T> *)
  {
    types::dynamic_tuple<T> res;
    res.data->resize(n);
    return res;
  }
  template <class T, size_t N>
  types::array_tuple<T, N> init_permut_from(size_t n, types::array_tuple<T, N> *)
  {
    assert(N == n && "consistent init");
    return {};
  }

  template <class T, class H>
  H permutations_iterator<T, H>::operator*() const
  {
    H res = init_permut_from(_size, (H *)nullptr);
    for (size_t i = 0; i < _size; i++)
      res[i] = pool[curr_permut[i]]; // Ok because types::dynamic_tuple is
                                     // indeed a vector
    return res;
  }

  template <class T, class I>
  types::dynamic_tuple<T> init_permut_from(I begin, I end, types::dynamic_tuple<T> *)
  {
    return {begin, end};
  }
  template <class T, size_t N, class I>
  types::array_tuple<T, N> init_permut_from(I begin, I end, types::array_tuple<T, N> *)
  {
    types::array_tuple<T, N> res;
    std::copy(begin, end, res.begin());
    return res;
  }

  template <class T, class H>
  permutations_iterator<T, H> &permutations_iterator<T, H>::operator++()
  {
    if (_size != pool.size()) {
      // Slow path, the iterator is a "view" of a prefix smaller
      // than the the pool size
      // FIXME a better implementation would be to avoid
      // std::next_permutation, but only in the slow path
      H prev_permut =
          init_permut_from(curr_permut.begin(), curr_permut.begin() + _size, (H *)nullptr);
      while ((end = std::next_permutation(curr_permut.begin(), curr_permut.end()))) {
        // Check if the prefix of the new permutation is
        // different of the previous one
        H new_permut =
            init_permut_from(curr_permut.begin(), curr_permut.begin() + _size, (H *)nullptr);
        if (!(prev_permut == new_permut))
          break;
      }
    } else
      end = std::next_permutation(curr_permut.begin(), curr_permut.end());
    return *this;
  }

  template <class T, class H>
  bool permutations_iterator<T, H>::operator!=(permutations_iterator<T, H> const &other) const
  {
    return !(*this == other);
  }

  template <class T, class H>
  bool permutations_iterator<T, H>::operator==(permutations_iterator<T, H> const &other) const
  {
    if (other.end != end)
      return false;
    return std::equal(curr_permut.begin(), curr_permut.end(), other.curr_permut.begin());
  }

  template <class T, class H>
  bool permutations_iterator<T, H>::operator<(permutations_iterator<T, H> const &other) const
  {
    if (end != other.end)
      return end > other.end;
    for (long i = 0; i < pool.size(); i++)
      if (other.curr_permut[i] < curr_permut[i])
        return false;
      else if (other.curr_permut[i] > curr_permut[i])
        return true;
    return false;
  }

  template <class T, class H>
  _permutations<T, H>::_permutations()
  {
  }

  template <class T, class H>
  _permutations<T, H>::_permutations(T iter, long elts)
      : iterator({iter.begin(), iter.end()}, elts, true)
  {
  }

  template <class T, class H>
  typename _permutations<T, H>::iterator const &_permutations<T, H>::begin() const
  {
    return *this;
  }

  template <class T, class H>
  typename _permutations<T, H>::iterator _permutations<T, H>::begin()
  {
    return *this;
  }

  template <class T, class H>
  typename _permutations<T, H>::iterator _permutations<T, H>::end() const
  {
    return iterator(iterator::pool, iterator::_size, false);
  }

  template <typename T0>
  _permutations<T0, types::dynamic_tuple<typename T0::value_type>> permutations(T0 iter,
                                                                                long num_elts)
  {
    return {iter, num_elts};
  }

  template <typename T0>
  _permutations<T0, types::dynamic_tuple<typename T0::value_type>> permutations(T0 iter)
  {
    return {iter, std::distance(iter.begin(), iter.end())};
  }

  template <typename T, long N>
  _permutations<T, types::array_tuple<typename T::value_type, (size_t)N>>
  permutations(T iter, std::integral_constant<long, N>)
  {
    return {iter, N};
  }
} // namespace itertools
PYTHONIC_NS_END

#endif
