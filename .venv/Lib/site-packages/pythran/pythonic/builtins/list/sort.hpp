#ifndef PYTHONIC_BUILTIN_LIST_SORT_HPP
#define PYTHONIC_BUILTIN_LIST_SORT_HPP

#include "pythonic/include/builtins/list/sort.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/types/NoneType.hpp"
#include "pythonic/types/list.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/pdqsort.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace list
  {

    template <class T>
    types::none_type sort(types::list<T> &seq)
    {
      pdqsort(seq.begin(), seq.end());
      return builtins::None;
    }

    template <class T, class K>
    types::none_type sort(types::list<T> &seq, K key)
    {
      pdqsort(seq.begin(), seq.end(),
              [&key](T const &self, T const &other) { return key(self) < key(other); });
      return builtins::None;
    }
  } // namespace list
} // namespace builtins
PYTHONIC_NS_END
#endif
