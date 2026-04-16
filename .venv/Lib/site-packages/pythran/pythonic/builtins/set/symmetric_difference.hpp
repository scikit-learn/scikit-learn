#ifndef PYTHONIC_BUILTIN_SET_SYMMETRICDIFFERENCE_HPP
#define PYTHONIC_BUILTIN_SET_SYMMETRICDIFFERENCE_HPP

#include "pythonic/include/builtins/set/symmetric_difference.hpp"

#include "pythonic/types/set.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace set
  {

    template <typename T, typename U>
    typename __combined<types::set<T>, U>::type symmetric_difference(types::set<T> const &set,
                                                                     U const &other)
    {
      return set.symmetric_difference(other);
    }

    /* No rvalue overload possible because of return type modification.:
     * >>> a = set([1, 2, 3])
     * >>> b = set([2., 3., 4.])
     * >>> a.symmetric_difference(b)
     * set([1.0, 4.0])
     */
    template <typename U>
    typename __combined<types::empty_set, U>::type symmetric_difference(types::empty_set const &set,
                                                                        U const &other)
    {
      return other;
    }
  } // namespace set
} // namespace builtins
PYTHONIC_NS_END
#endif
