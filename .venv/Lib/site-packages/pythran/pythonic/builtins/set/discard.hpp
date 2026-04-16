#ifndef PYTHONIC_BUILTIN_SET_DISCARD_HPP
#define PYTHONIC_BUILTIN_SET_DISCARD_HPP

#include "pythonic/include/builtins/set/discard.hpp"

#include "pythonic/types/set.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {
    template <class T, class U>
    void discard(types::set<T> &set, U const &elem)
    {
      set.discard(elem);
    }

    template <class T, class U>
    void discard(types::set<T> &&set, U const &elem)
    {
      // nothing to be done for lvalue
    }

    template <class U>
    void discard(types::empty_set const &set, U const &elem)
    {
      // nothing to remove in an empty_set
    }
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
