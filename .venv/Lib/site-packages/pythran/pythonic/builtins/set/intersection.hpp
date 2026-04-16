#ifndef PYTHONIC_BUILTIN_SET_INTERSECTION_HPP
#define PYTHONIC_BUILTIN_SET_INTERSECTION_HPP

#include "pythonic/include/builtins/set/intersection.hpp"

#include "pythonic/types/set.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {

    template <typename T, typename... Types>
    typename __combined<types::set<T>, Types...>::type intersection(types::set<T> const &set,
                                                                    Types const &...others)
    {
      return set.intersection(others...);
    }

    /* No rvalue overload possible because of return type modification.:
     * >>> a = set([1,2,3])
     * >>> b = set([1., 2., 3.])
     * >>> a.intersection(b)
     * set([1.0, 2.0, 3.0])
     */
    template <typename... Types>
    types::empty_set intersection(types::empty_set const &set, Types const &...others)
    {
      return types::empty_set();
    }
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
