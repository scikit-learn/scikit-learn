#ifndef PYTHONIC_INCLUDE_ITERTOOLS_COUNT_HPP
#define PYTHONIC_INCLUDE_ITERTOOLS_COUNT_HPP

#include "pythonic/include/types/combined.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <iterator>

PYTHONIC_NS_BEGIN

namespace itertools
{
  namespace details
  {
    template <class T>
    struct count_iterator : std::iterator<std::random_access_iterator_tag, T> {
      T value;
      T step;
      count_iterator() = default;
      count_iterator(T value, T step);
      T operator*() const;
      count_iterator &operator++();
      count_iterator &operator+=(long n);
      bool operator!=(count_iterator const &other) const;
      bool operator==(count_iterator const &other) const;
      bool operator<(count_iterator const &other) const;
      long operator-(count_iterator const &other) const;
    };

    template <class T>
    struct count : count_iterator<T> {
      using value_type = T;
      using iterator = count_iterator<T>;

      count() = default;
      count(T value, T step);
      iterator &begin();
      iterator const &begin() const;
      iterator end() const;
    };
  } // namespace details

  template <typename T0, typename T1 = T0>
  details::count<typename __combined<T0, T1>::type> count(T0 start, T1 step = 1);

  details::count<long> count();

  DEFINE_FUNCTOR(pythonic::itertools, count);
} // namespace itertools
PYTHONIC_NS_END

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class E, class T>
struct __combined<E, pythonic::itertools::details::count<T>> {
  using type = typename __combined<
      E, container<typename pythonic::itertools::details::count<T>::value_type>>::type;
};

/* } */

#endif
