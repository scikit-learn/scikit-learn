#ifndef PYTHONIC_INCLUDE_ITERTOOLS_ISLICE_HPP
#define PYTHONIC_INCLUDE_ITERTOOLS_ISLICE_HPP

#include "pythonic/include/builtins/range.hpp"
#include "pythonic/include/itertools/common.hpp"
#include "pythonic/include/utils/functor.hpp"
#include <iterator>

PYTHONIC_NS_BEGIN

namespace itertools
{
  template <typename Iterable>
  struct islice_iterator
      : std::iterator<typename Iterable::iterator::iterator_category,
                      typename std::iterator_traits<typename Iterable::iterator>::value_type> {
    std::remove_reference_t<std::remove_cv_t<Iterable>> iterable_ref;
    typename std::remove_reference_t<std::remove_cv_t<Iterable>>::iterator iterable;

    builtins::range xr_ref;
    builtins::range_iterator state;
    builtins::range_iterator::value_type prev;

    islice_iterator();
    islice_iterator(Iterable const &iterable, builtins::range const &xr);
    islice_iterator(npos const &n, Iterable const &iterable, builtins::range const &xr);

    typename Iterable::value_type operator*() const;
    islice_iterator &operator++();
    bool operator==(islice_iterator const &other) const;
    bool operator!=(islice_iterator const &other) const;
    bool operator<(islice_iterator const &other) const;
    int operator-(islice_iterator const &other) const;
  };

  template <typename Iterable>
  struct _islice : islice_iterator<Iterable> {

    using iterator = islice_iterator<Iterable>;
    using value_type = typename Iterable::value_type;

    iterator end_iter;

    _islice();
    _islice(Iterable const &iterable, builtins::range const &xr);

    iterator &begin();
    iterator const &begin() const;
    iterator end() const;
  };

  template <typename Iterable>
  _islice<std::remove_cv_t<std::remove_reference_t<Iterable>>>
  islice(Iterable &&iterable, long start, long stop, long step = 1);

  template <typename Iterable>
  _islice<std::remove_cv_t<std::remove_reference_t<Iterable>>> islice(Iterable &&iterable,
                                                                      long stop);

  DEFINE_FUNCTOR(pythonic::itertools, islice);
} // namespace itertools
PYTHONIC_NS_END

/* type inference stuff  {*/
#include "pythonic/include/types/combined.hpp"

template <class E, class T>
struct __combined<E, pythonic::itertools::_islice<T>> {
  using type =
      typename __combined<E, container<typename pythonic::itertools::_islice<T>::value_type>>::type;
};

/* } */

#endif
