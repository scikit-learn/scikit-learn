#ifndef PYTHONIC_ITERTOOLS_ISLICE_HPP
#define PYTHONIC_ITERTOOLS_ISLICE_HPP

#include "pythonic/builtins/range.hpp"
#include "pythonic/include/itertools/islice.hpp"
#include "pythonic/itertools/common.hpp"
#include "pythonic/utils/functor.hpp"
#include <iterator>

PYTHONIC_NS_BEGIN

namespace itertools
{
  template <typename Iterable>
  islice_iterator<Iterable>::islice_iterator()
  {
  }

  template <typename Iterable>
  islice_iterator<Iterable>::islice_iterator(Iterable const &iterable, builtins::range const &xr)
      : iterable_ref(iterable), iterable(iterable_ref.begin()), xr_ref(xr), state(xr_ref.begin()),
        prev(*state)
  {
    std::advance(this->iterable, *state);
  }

  template <typename Iterable>
  islice_iterator<Iterable>::islice_iterator(npos const &n, Iterable const &iterable,
                                             builtins::range const &xr)
      : iterable_ref(iterable), iterable(iterable_ref.begin()), xr_ref(xr), state(xr_ref.end()),
        prev(0)
  {
  }

  template <typename Iterable>
  typename Iterable::value_type islice_iterator<Iterable>::operator*() const
  {
    return *iterable;
  }

  template <typename Iterable>
  islice_iterator<Iterable> &islice_iterator<Iterable>::operator++()
  {
    ++state;
    std::advance(this->iterable, *state - prev);
    prev = *state;
    return *this;
  }

  template <typename Iterable>
  bool islice_iterator<Iterable>::operator==(islice_iterator<Iterable> const &other) const
  {
    return (state == other.state);
  }

  template <typename Iterable>
  bool islice_iterator<Iterable>::operator!=(islice_iterator<Iterable> const &other) const
  {
    return state != other.state;
  }

  template <typename Iterable>
  bool islice_iterator<Iterable>::operator<(islice_iterator<Iterable> const &other) const
  {
    return state != other.state;
  }

  template <typename Iterable>
  int islice_iterator<Iterable>::operator-(islice_iterator<Iterable> const &other) const
  {
    return state - other.state;
  }

  template <typename Iterable>
  _islice<Iterable>::_islice()
  {
  }

  template <typename Iterable>
  _islice<Iterable>::_islice(Iterable const &iterable, builtins::range const &xr)
      : iterator(iterable, xr), end_iter(npos(), iterable, xr)
  {
  }

  template <typename Iterable>
  typename _islice<Iterable>::iterator &_islice<Iterable>::begin()
  {
    return *this;
  }

  template <typename Iterable>
  typename _islice<Iterable>::iterator const &_islice<Iterable>::begin() const
  {
    return *this;
  }

  template <typename Iterable>
  typename _islice<Iterable>::iterator _islice<Iterable>::end() const
  {
    return end_iter;
  }

  template <typename Iterable>
  _islice<std::remove_cv_t<std::remove_reference_t<Iterable>>>
  islice(Iterable &&iterable, long start, long stop, long step)
  {
    return {iterable, builtins::range(start, stop, step)};
  }

  template <typename Iterable>
  _islice<std::remove_cv_t<std::remove_reference_t<Iterable>>> islice(Iterable &&iterable,
                                                                      long stop)
  {
    return {iterable, builtins::range(0, stop, 1)};
  }
} // namespace itertools
PYTHONIC_NS_END

#endif
