#ifndef PYTHONIC_ITERTOOLS_REPEAT_HPP
#define PYTHONIC_ITERTOOLS_REPEAT_HPP

#include "pythonic/include/itertools/repeat.hpp"
#include "pythonic/types/list.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace itertools
{

  template <class T, bool Endless>
  repeat_iterator<T, Endless>::repeat_iterator(T value, long count) : value_(value), count_(count)
  {
  }

  template <class T, bool Endless>
  repeat_iterator<T, Endless> &repeat_iterator<T, Endless>::operator++()
  {
    ++count_;
    return *this;
  }

  template <class T, bool Endless>
  T repeat_iterator<T, Endless>::operator*()
  {
    return value_;
  }

  template <class T, bool Endless>
  bool repeat_iterator<T, Endless>::operator!=(repeat_iterator<T, Endless> const &other) const
  {
    return Endless || count_ != other.count_;
  }
  template <class T, bool Endless>
  bool repeat_iterator<T, Endless>::operator==(repeat_iterator<T, Endless> const &other) const
  {
    return !Endless && count_ == other.count_;
  }
  template <class T, bool Endless>
  bool repeat_iterator<T, Endless>::operator<(repeat_iterator<T, Endless> const &other) const
  {
    return !Endless && count_ < other.count_;
  }

  template <class T, bool Endless>
  _repeat<T, Endless>::_repeat(T value, long count) : repeat_iterator<T, Endless>(value, count)
  {
  }

  template <class T, bool Endless>
  typename _repeat<T, Endless>::iterator _repeat<T, Endless>::begin() const
  {
    return {_repeat<T, Endless>::iterator::value_, 0};
  }
  template <class T, bool Endless>
  typename _repeat<T, Endless>::iterator _repeat<T, Endless>::end() const
  {
    return *this;
  }

  template <typename T>
  _repeat<T, false> repeat(T value, long count)
  {
    return {value, count};
  }

  template <typename T>
  _repeat<T, true> repeat(T value)
  {
    return {value, -1};
  }
} // namespace itertools
PYTHONIC_NS_END

#endif
