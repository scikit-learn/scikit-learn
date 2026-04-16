#ifndef PYTHONIC_INCLUDE_BUILTIN_RANGE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_RANGE_HPP

#include "pythonic/include/types/list.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  namespace
  {
    struct range_iterator : std::iterator<std::random_access_iterator_tag, long, ptrdiff_t, long *,
                                          long /*no ref here*/> {
      long value_;
      long step_;

      range_iterator() = default;
      range_iterator(long v, long s);
      long operator*() const;
      range_iterator &operator++();
      range_iterator &operator--();
      range_iterator operator++(int);
      range_iterator operator--(int);
      range_iterator &operator+=(long n);
      range_iterator &operator-=(long n);
      bool operator!=(range_iterator const &other) const;
      bool operator==(range_iterator const &other) const;
      bool operator<(range_iterator const &other) const;
      bool operator<=(range_iterator const &other) const;
      long operator-(range_iterator const &other) const;
    };
  } // namespace

  struct range {
    using value_type = long;
    using iterator = range_iterator;
    using const_iterator = range_iterator;
    using reverse_iterator = range_iterator;
    using const_reverse_iterator = range_iterator;
    using dtype = long;
    static constexpr long value = 1;

    long begin_;
    long end_;
    long step_;

    range() = default;
    range(long b, long e, long s = 1);
    range(long e);
    iterator begin() const;
    iterator end() const;
    reverse_iterator rbegin() const;
    reverse_iterator rend() const;

    long size() const;
    long operator[](long i) const;
  };

  DEFINE_FUNCTOR(pythonic::builtins, range);
} // namespace builtins
PYTHONIC_NS_END

namespace std
{
  template <size_t I>
  long get(pythonic::builtins::range const &);

  template <size_t I>
  struct tuple_element<I, pythonic::builtins::range> {
    typedef long type;
  };
} // namespace std

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class E>
struct __combined<E, pythonic::builtins::range> {
  using type = typename __combined<E, container<long>>::type;
};

#endif
