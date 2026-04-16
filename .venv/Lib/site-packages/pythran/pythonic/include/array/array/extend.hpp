#ifndef PYTHONIC_INCLUDE_ARRAY_ARRAY_EXTEND_HPP
#define PYTHONIC_INCLUDE_ARRAY_ARRAY_EXTEND_HPP

#include "pythonic/include/types/array.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace array
{

  namespace array
  {

    template <class T, class S>
    types::none_type extend(types::array<T> &a, S &&elts);

    DEFINE_FUNCTOR(pythonic::array::array, extend);
  } // namespace array
} // namespace array
PYTHONIC_NS_END
#endif
