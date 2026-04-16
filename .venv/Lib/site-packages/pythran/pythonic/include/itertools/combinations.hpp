#ifndef PYTHONIC_INCLUDE_ITERTOOLS_COMBINATIONS_HPP
#define PYTHONIC_INCLUDE_ITERTOOLS_COMBINATIONS_HPP

#include "pythonic/include/types/dynamic_tuple.hpp"
#include "pythonic/include/utils/allocate.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <iterator>
#include <vector>

PYTHONIC_NS_BEGIN

namespace itertools
{
  namespace details
  {
    template <class T>
    struct combination_iterator
        : std::iterator<std::forward_iterator_tag, types::dynamic_tuple<typename T::value_type>,
                        ptrdiff_t, types::dynamic_tuple<typename T::value_type> *,
                        types::dynamic_tuple<typename T::value_type> /*no ref*/
                        > {
      std::vector<typename T::value_type, utils::allocator<typename T::value_type>> pool;
      std::vector<long, utils::allocator<long>> indices;
      long r;
      bool stopped;
      std::vector<typename T::value_type, utils::allocator<typename T::value_type>> result;

      combination_iterator() = default;
      combination_iterator(bool);

      template <class Iter>
      combination_iterator(Iter &&pool, long r);

      types::dynamic_tuple<typename T::value_type> operator*() const;
      combination_iterator &operator++();
      bool operator!=(combination_iterator const &other) const;
      bool operator==(combination_iterator const &other) const;
      bool operator<(combination_iterator const &other) const;
    };

    template <class T>
    struct combination : combination_iterator<T> {
      using iterator = combination_iterator<T>;
      using value_type = typename iterator::value_type;

      long num_elts;

      combination() = default;

      template <class Iter>
      combination(Iter &&iter, long elts);
      iterator const &begin() const;
      iterator begin();
      iterator end() const;
    };
  } // namespace details

  template <typename T0>
  details::combination<std::remove_cv_t<std::remove_reference_t<T0>>> combinations(T0 &&iter,
                                                                                   long num_elts);

  DEFINE_FUNCTOR(pythonic::itertools, combinations);
} // namespace itertools
PYTHONIC_NS_END

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class E, class T>
struct __combined<E, pythonic::itertools::details::combination<T>> {
  using type = typename __combined<
      E, container<typename pythonic::itertools::details::combination<T>::value_type>>::type;
};

/* } */

#endif
