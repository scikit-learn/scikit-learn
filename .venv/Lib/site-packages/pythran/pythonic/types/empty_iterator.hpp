#ifndef PYTHONIC_TYPES_EMPTY_ITERATOR_HPP
#define PYTHONIC_TYPES_EMPTY_ITERATOR_HPP

#include "pythonic/include/types/empty_iterator.hpp"

PYTHONIC_NS_BEGIN

namespace types
{

  empty_iterator::empty_iterator()
  {
  }

  empty_iterator::empty_iterator(empty_iterator const &)
  {
  }

  inline bool empty_iterator::operator==(empty_iterator const &) const
  {
    return true;
  }

  inline bool empty_iterator::operator!=(empty_iterator const &) const
  {
    return false;
  }

  inline bool empty_iterator::operator<(empty_iterator const &) const
  {
    return false;
  }

  inline empty_iterator &empty_iterator::operator++()
  {
    return *this;
  }

  inline empty_iterator &empty_iterator::operator++(int)
  {
    return *this;
  }

  inline double empty_iterator::operator*() const
  {
    return {};
  }

  inline void empty_iterator::operator->() const
  {
    return;
  }
} // namespace types
PYTHONIC_NS_END

#endif
