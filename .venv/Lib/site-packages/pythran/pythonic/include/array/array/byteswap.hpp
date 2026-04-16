#ifndef PYTHONIC_INCLUDE_ARRAY_ARRAY_BYTESWAP_HPP
#define PYTHONIC_INCLUDE_ARRAY_ARRAY_BYTESWAP_HPP

#include "pythonic/include/types/array.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace array
{

  namespace array
  {

    template <class T>
    types::none_type byteswap(types::array<T> &seq);

    template <class T>
    types::none_type byteswap(types::array<T> &&seq)
    {
      return {};
    }

    DEFINE_FUNCTOR(pythonic::array::array, byteswap);
  } // namespace array
} // namespace array
PYTHONIC_NS_END
#endif
