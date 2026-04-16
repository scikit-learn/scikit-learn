#ifndef PYTHONIC_ARRAY_ARRAY_BUFFER_INFO_HPP
#define PYTHONIC_ARRAY_ARRAY_BUFFER_INFO_HPP

#include "pythonic/include/types/array.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace array
{

  namespace array
  {

    template <class T>
    std::tuple<long, long> buffer_info(types::array<T> const &seq)
    {
      return std::make_tuple(seq.id(), seq.size());
    }

  } // namespace array
} // namespace array
PYTHONIC_NS_END
#endif
