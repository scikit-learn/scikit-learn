#ifndef PYTHONIC_INCLUDE_BUILTIN_XRANGE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_XRANGE_HPP

#include "pythonic/include/utils/functor.hpp"
#include <iterator>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace
  {
    struct xrange_iterator : std::iterator<std::random_access_iterator_tag, long, ptrdiff_t, long *,
                                           long /*no ref here*/> {
      long value_;
      long step_;

      xrange_iterator() = default;
      xrange_iterator(long v, long s);
      long operator*() const;
      xrange_iterator &operator++();
      xrange_iterator operator++(int);
      xrange_iterator &operator+=(long n);
      bool operator!=(xrange_iterator const &other) const;
      bool operator==(xrange_iterator const &other) const;
      bool operator<(xrange_iterator const &other) const;
      long operator-(xrange_iterator const &other) const;
    };
  } // namespace

  struct xrange {
    using value_type = long;
    using iterator = xrange_iterator;
    using const_iterator = xrange_iterator;
    using reverse_iterator = xrange_iterator;
    using const_reverse_iterator = xrange_iterator;

    long begin_;
    long end_;
    long step_;

    xrange() = default;
    xrange(long b, long e, long s = 1);
    xrange(long e);
    iterator begin() const;
    iterator end() const;
    reverse_iterator rbegin() const;
    reverse_iterator rend() const;
  };

  DEFINE_FUNCTOR(pythonic::builtins, xrange);
} // namespace builtins
PYTHONIC_NS_END

#endif
