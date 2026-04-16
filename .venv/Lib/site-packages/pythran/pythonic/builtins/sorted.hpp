#ifndef PYTHONIC_BUILTIN_SORTED_HPP
#define PYTHONIC_BUILTIN_SORTED_HPP

#include "pythonic/include/builtins/sorted.hpp"

#include "pythonic/types/list.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/pdqsort.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class Iterable>
  types::list<std::remove_cv_t<
      typename std::iterator_traits<typename std::decay_t<Iterable>::iterator>::value_type>>
  sorted(Iterable &&seq)
  {
    types::list<std::remove_cv_t<
        typename std::iterator_traits<typename std::decay_t<Iterable>::iterator>::value_type>>
        out(seq.begin(), seq.end());
    pdqsort(out.begin(), out.end());
    return out;
  }

  template <class Iterable, class Key>
  types::list<std::remove_cv_t<
      typename std::iterator_traits<typename std::decay_t<Iterable>::iterator>::value_type>>
  sorted(Iterable &&seq, types::kwonly, Key const &key, bool reverse)
  {
    using value_type = std::remove_cv_t<
        typename std::iterator_traits<typename std::decay_t<Iterable>::iterator>::value_type>;
    types::list<value_type> out(seq.begin(), seq.end());
    if (reverse)
      pdqsort(out.begin(), out.end(), [&key](value_type const &self, value_type const &other) {
        return key(self) > key(other);
      });
    else
      pdqsort(out.begin(), out.end(), [&key](value_type const &self, value_type const &other) {
        return key(self) < key(other);
      });
    return out;
  }

  template <class Iterable>
  types::list<std::remove_cv_t<
      typename std::iterator_traits<typename std::decay_t<Iterable>::iterator>::value_type>>
  sorted(Iterable &&seq, types::kwonly, types::none_type const &key, bool reverse)
  {
    using value_type = std::remove_cv_t<
        typename std::iterator_traits<typename std::decay_t<Iterable>::iterator>::value_type>;
    types::list<value_type> out(seq.begin(), seq.end());
    if (reverse)
      pdqsort(out.begin(), out.end(),
              [](value_type const &self, value_type const &other) { return self > other; });
    else
      pdqsort(out.begin(), out.end());
    return out;
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
