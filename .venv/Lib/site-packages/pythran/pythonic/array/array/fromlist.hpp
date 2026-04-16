#ifndef PYTHONIC_ARRAY_ARRAY_FROMLIST_HPP
#define PYTHONIC_ARRAY_ARRAY_FROMLIST_HPP

#include "pythonic/include/types/array.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace array
{

  namespace array
  {

    template <class T, class E>
    types::none_type fromlist(types::array<T> &seq, E &&elts)
    {
      seq += std::forward<E>(elts);
      return {};
    }

  } // namespace array
} // namespace array
PYTHONIC_NS_END
#endif
