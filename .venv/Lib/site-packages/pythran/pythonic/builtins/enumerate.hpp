#ifndef PYTHONIC_BUILTIN_ENUMERATE_HPP
#define PYTHONIC_BUILTIN_ENUMERATE_HPP

#include "pythonic/include/builtins/enumerate.hpp"

#include "pythonic/utils/functor.hpp"

#include <tuple>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace details
  {
    /// enumerate_iterator implementation

    template <class Iterator>
    enumerate_iterator<Iterator>::enumerate_iterator()
    {
    }

    template <class Iterator>
    enumerate_iterator<Iterator>::enumerate_iterator(Iterator const &iter, long first)
        : value(first), iter(iter)
    {
    }

    template <class Iterator>
    enumerate_iterator<Iterator> &enumerate_iterator<Iterator>::operator+=(long n)
    {
      value += n, iter += n;
      return *this;
    }

    // Comparison operators can't use value as end() doesn't have a valid
    // value content
    // du to the lake of size information for generator
    // TODO : We could handle case with && without size if there is a
    // performances benefits
    template <class Iterator>
    bool enumerate_iterator<Iterator>::operator!=(enumerate_iterator<Iterator> const &other) const
    {
      return !(*this == other);
    }

    template <class Iterator>
    bool enumerate_iterator<Iterator>::operator<(enumerate_iterator const &other) const
    {
      return iter < other.iter;
    }

    template <class Iterator>
    bool enumerate_iterator<Iterator>::operator==(enumerate_iterator<Iterator> const &other) const
    {
      return iter == other.iter;
    }

    template <class Iterator>
    long enumerate_iterator<Iterator>::operator-(enumerate_iterator<Iterator> const &other) const
    {
      return iter - other.iter;
    }

    /// details::enumerate implementation
    template <class Iterable>
    enumerate<Iterable>::enumerate()
    {
    }

    template <class Iterable>
    enumerate<Iterable>::enumerate(Iterable seq, long first)
        : Iterable(seq), iterator(Iterable::begin(), first), end_iter(Iterable::end(), -1)
    {
    }

    template <class Iterable>
    typename enumerate<Iterable>::iterator &enumerate<Iterable>::begin()
    {
      return *this;
    }

    template <class Iterable>
    typename enumerate<Iterable>::iterator const &enumerate<Iterable>::begin() const
    {
      return *this;
    }

    template <class Iterable>
    typename enumerate<Iterable>::iterator enumerate<Iterable>::end() const
    {
      return end_iter;
    }
  } // namespace details

  /// enumerate implementation

  template <class Iterable>
  details::enumerate<std::remove_cv_t<std::remove_reference_t<Iterable>>> enumerate(Iterable &&seq,
                                                                                    long first)
  {
    return {std::forward<Iterable>(seq), first};
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
