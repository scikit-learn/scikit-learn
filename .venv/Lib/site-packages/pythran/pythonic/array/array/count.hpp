#ifndef PYTHONIC_ARRAY_ARRAY_COUNT_HPP
#define PYTHONIC_ARRAY_ARRAY_COUNT_HPP

#include "pythonic/include/types/array.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace array
{

  namespace array
  {

    template <class T>
    long count(types::array<T> const &seq)
    {
      return std::count(seq.begin(), seq.end());
    }

  } // namespace array
} // namespace array
PYTHONIC_NS_END
#endif
