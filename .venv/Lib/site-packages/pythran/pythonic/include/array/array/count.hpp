#ifndef PYTHONIC_INCLUDE_ARRAY_ARRAY_COUNT_HPP
#define PYTHONIC_INCLUDE_ARRAY_ARRAY_COUNT_HPP

#include "pythonic/include/types/array.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace array
{

  namespace array
  {

    template <class T>
    long count(types::array<T> const &seq);

    DEFINE_FUNCTOR(pythonic::array::array, count);
  } // namespace array
} // namespace array
PYTHONIC_NS_END
#endif
