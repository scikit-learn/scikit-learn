#ifndef PYTHONIC_INCLUDE_ITERTOOLS_PERMUTATIONS_HPP
#define PYTHONIC_INCLUDE_ITERTOOLS_PERMUTATIONS_HPP

#include "pythonic/include/types/dynamic_tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

#include "pythonic/utils/allocate.hpp"

#include <iterator>
#include <vector>

PYTHONIC_NS_BEGIN

namespace itertools
{
  /** Permutation iterator
   *
   *  It wraps a vector && provide an iteration over every possible
   *  permutation of the vector. The permutations are represented as
   *dynamic_tuple
   *  of elements.
   *
   *  The following iterator:
   *
   *  permutations_iterator([0, 1, 2])
   *
   *  yields the following suite:
   *
   *  [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
   *
   */
  template <class T, class H>
  struct permutations_iterator
      : std::iterator<std::forward_iterator_tag, H, ptrdiff_t, H *, H /* no ref*/
                      > {

    // Vector of inputs, contains elements to permute
    using pool_type = std::vector<typename T::value_type, utils::allocator<typename T::value_type>>;
    pool_type pool;

    // The current permutation as a dynamic_tuple of index in the pool
    // Internally it always has the same size as the pool, even if the
    // external view is limited
    std::vector<long, utils::allocator<long>> curr_permut;

    // Size of the "visible" permutation
    size_t _size;
    bool end; // sentinel marker

    permutations_iterator();
    permutations_iterator(pool_type const &iter, size_t num_elts, bool end);

    /** Build the permutation visible from the "outside" */
    H operator*() const;

    /*  Generate next permutation
     *
     *  If the size of the permutation is smaller than the size of the
     *  pool, we may have to iterate multiple times
     */
    permutations_iterator &operator++();
    bool operator!=(permutations_iterator const &other) const;
    bool operator==(permutations_iterator const &other) const;
    bool operator<(permutations_iterator const &other) const;
  };

  template <class T, class H>
  // FIXME document why this inheritance???
  struct _permutations : permutations_iterator<T, H> {
    using iterator = permutations_iterator<T, H>;
    using value_type = typename iterator::value_type;

    _permutations();
    _permutations(T iter, long elts);

    iterator const &begin() const;
    iterator begin();
    iterator end() const;
  };

  template <typename T0>
  _permutations<T0, types::dynamic_tuple<typename T0::value_type>> permutations(T0 iter,
                                                                                long num_elts);

  template <typename T0>
  _permutations<T0, types::dynamic_tuple<typename T0::value_type>> permutations(T0 iter);

  template <typename T0, long N0>
  _permutations<T0, types::array_tuple<typename T0::value_type, (size_t)N0>>
  permutations(T0 iter, std::integral_constant<long, N0>);

  DEFINE_FUNCTOR(pythonic::itertools, permutations);
} // namespace itertools
PYTHONIC_NS_END

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class E, class T, class H>
struct __combined<E, pythonic::itertools::_permutations<T, H>> {
  using type = typename __combined<
      E, container<typename pythonic::itertools::_permutations<T, H>::value_type>>::type;
};

/* } */

#endif
