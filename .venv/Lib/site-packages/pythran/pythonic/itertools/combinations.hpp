#ifndef PYTHONIC_ITERTOOLS_COMBINATIONS_HPP
#define PYTHONIC_ITERTOOLS_COMBINATIONS_HPP

#include "pythonic/include/itertools/combinations.hpp"

#include "pythonic/types/dynamic_tuple.hpp"
#include "pythonic/utils/functor.hpp"

#include <numeric>

PYTHONIC_NS_BEGIN

namespace itertools
{
  namespace details
  {
    template <class T>
    template <class Iter>
    combination_iterator<T>::combination_iterator(Iter &&pool, long r)
        : pool(pool.begin(), pool.end()), indices(r), r(r), stopped(r > long(this->pool.size()))
    {
      assert(r >= 0 && "r must be non-negative");
      if (!stopped) {
        std::iota(indices.begin(), indices.end(), 0);
        result.insert(result.end(), this->pool.begin(), this->pool.begin() + r);
      }
    }

    template <class T>
    combination_iterator<T>::combination_iterator(bool) : stopped(true)
    {
    }

    template <class T>
    types::dynamic_tuple<typename T::value_type> combination_iterator<T>::operator*() const
    {
      assert(!stopped && "! stopped");
      return {result.begin(), result.end()};
    }

    template <class T>
    combination_iterator<T> &combination_iterator<T>::operator++()
    {
      /* Scan indices right-to-left until finding one that is !
         at its maximum (i + n - r). */
      long i, n = pool.size();
      for (i = r - 1; i >= 0 && indices[i] == i + n - r; i--)
        ;

      /* If i is negative, then the indices are all at
         their maximum value && we're done. */
      if (i < 0)
        stopped = true;
      else {
        /* Increment the current index which we know is ! at its
           maximum.  Then move back to the right setting each index
           to its lowest possible value (one higher than the index
           to its left -- this maintains the sort order invariant). */
        indices[i]++;
        for (long j = i + 1; j < r; j++)
          indices[j] = indices[j - 1] + 1;

        /* Update the result tuple for the new indices
           starting with i, the leftmost index that changed */
        for (; i < r; i++) {
          result[i] = pool[indices[i]];
        }
      }
      return *this;
    }

    template <class T>
    bool combination_iterator<T>::operator!=(combination_iterator const &other) const
    {
      assert(stopped || other.stopped);
      return !(*this == other);
    }

    template <class T>
    bool combination_iterator<T>::operator==(combination_iterator const &other) const
    {
      assert(stopped || other.stopped);
      return other.stopped == stopped;
    }

    template <class T>
    bool combination_iterator<T>::operator<(combination_iterator const &other) const
    {
      return stopped != other.stopped;
    }

    template <class T>
    template <class Iter>
    combination<T>::combination(Iter &&iter, long elts)
        : iterator(std::forward<Iter>(iter), elts), num_elts(elts)
    {
    }

    template <class T>
    typename combination<T>::iterator const &combination<T>::begin() const
    {
      return *this;
    }

    template <class T>
    typename combination<T>::iterator combination<T>::begin()
    {
      return *this;
    }

    template <class T>
    typename combination<T>::iterator combination<T>::end() const
    {
      return {true};
    }
  } // namespace details

  template <typename T0>
  details::combination<std::remove_cv_t<std::remove_reference_t<T0>>> combinations(T0 &&iter,
                                                                                   long num_elts)
  {
    return {std::forward<T0>(iter), num_elts};
  }
} // namespace itertools
PYTHONIC_NS_END

#endif
