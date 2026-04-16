#ifndef PYTHONIC_UTILS_ITERATOR_HPP
#define PYTHONIC_UTILS_ITERATOR_HPP

#include "pythonic/include/utils/iterator.hpp"

PYTHONIC_NS_BEGIN

namespace utils
{

  template <class T>
  comparable_iterator<T>::comparable_iterator() : T()
  {
  }

  template <class T>
  comparable_iterator<T>::comparable_iterator(T const &t) : T(t)
  {
  }

  template <class T>
  bool comparable_iterator<T>::operator<(comparable_iterator<T> other)
  {
    return (*this) != other;
  }

  template <class T>
  iterator_reminder<false, T>::iterator_reminder(T const &v) : values(v)
  {
  }

  template <class T>
  iterator_reminder<true, T>::iterator_reminder(T const &v) : values(v)
  {
  }

  template <class T, class... Others>
  iterator_reminder<true, T, Others...>::iterator_reminder(T const &v, Others const &...others)
      : values(v, others...)
  {
  }
} // namespace utils
PYTHONIC_NS_END

#endif
