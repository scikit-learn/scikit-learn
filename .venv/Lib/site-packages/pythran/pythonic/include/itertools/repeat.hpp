#ifndef PYTHONIC_INCLUDE_ITERTOOLS_REPEAT_HPP
#define PYTHONIC_INCLUDE_ITERTOOLS_REPEAT_HPP

#include "pythonic/include/utils/functor.hpp"

#include <iterator>

PYTHONIC_NS_BEGIN

namespace itertools
{
  template <class T, bool Endless>
  struct repeat_iterator : std::iterator<std::forward_iterator_tag, T, ptrdiff_t, T *, T /* no ref*/
                                         > {
    T value_;
    long count_;

    repeat_iterator(T value, long count);
    repeat_iterator &operator++();
    T operator*();
    bool operator!=(repeat_iterator const &other) const;
    bool operator==(repeat_iterator const &other) const;
    bool operator<(repeat_iterator const &other) const;
  };

  template <class T, bool Endless>
  struct _repeat : repeat_iterator<T, Endless> {
    using iterator = repeat_iterator<T, Endless>;
    using value_type = typename iterator::value_type;

    _repeat() = default;
    _repeat(T value, long count);

    iterator begin() const;
    iterator end() const;
  };

  template <typename T>
  _repeat<T, false> repeat(T value, long num_elts);

  template <typename T>
  _repeat<T, true> repeat(T iter);

  DEFINE_FUNCTOR(pythonic::itertools, repeat);
} // namespace itertools
PYTHONIC_NS_END

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class E, class T, bool C>
struct __combined<E, pythonic::itertools::_repeat<T, C>> {
  using type =
      typename __combined<E,
                          container<typename pythonic::itertools::_repeat<T, C>::value_type>>::type;
};

/* } */

#endif
