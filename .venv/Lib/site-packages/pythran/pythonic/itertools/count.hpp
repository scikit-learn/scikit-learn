#ifndef PYTHONIC_ITERTOOLS_COUNT_HPP
#define PYTHONIC_ITERTOOLS_COUNT_HPP

#include "pythonic/include/itertools/count.hpp"

#include "pythonic/utils/functor.hpp"
#include <limits>

PYTHONIC_NS_BEGIN

namespace itertools
{
  namespace details
  {
    template <class T>
    count_iterator<T>::count_iterator(T value, T step) : value(value), step(step)
    {
    }

    template <class T>
    T count_iterator<T>::operator*() const
    {
      return value;
    }

    template <class T>
    count_iterator<T> &count_iterator<T>::operator++()
    {
      value += step;
      return *this;
    }

    template <class T>
    count_iterator<T> &count_iterator<T>::operator+=(long n)
    {
      value += step * n;
      return *this;
    }

    template <class T>
    bool count_iterator<T>::operator!=(count_iterator const &other) const
    {
      return value != other.value;
    }

    template <class T>
    bool count_iterator<T>::operator==(count_iterator const &other) const
    {
      return value == other.value;
    }

    template <class T>
    bool count_iterator<T>::operator<(count_iterator const &other) const
    {
      return value < other.value;
    }

    template <class T>
    long count_iterator<T>::operator-(count_iterator const &other) const
    {
      return (value - other.value) / step;
    }

    template <class T>
    count<T>::count(T value, T step) : count_iterator<T>(value, step)
    {
    }

    template <class T>
    typename count<T>::iterator &count<T>::begin()
    {
      return *this;
    }

    template <class T>
    typename count<T>::iterator const &count<T>::begin() const
    {
      return *this;
    }

    template <class T>
    typename count<T>::iterator count<T>::end() const
    {
      return {std::numeric_limits<T>::max(), count_iterator<T>::step};
    }
  } // namespace details

  template <typename T0, typename T1>
  details::count<typename __combined<T0, T1>::type> count(T0 start, T1 step)
  {
    using return_t = typename __combined<T0, T1>::type;
    return {static_cast<return_t>(start), static_cast<return_t>(step)};
  }

  inline details::count<long> count()
  {
    return {0, 1};
  }
} // namespace itertools
PYTHONIC_NS_END

#endif
