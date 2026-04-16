#ifndef PYTHONIC_INCLUDE_ARRAY_ARRAY_FROMBYTES_HPP
#define PYTHONIC_INCLUDE_ARRAY_ARRAY_FROMBYTES_HPP

#include "pythonic/include/types/array.hpp"
#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace array
{

  namespace array
  {

    template <class T>
    types::none_type frombytes(types::array<T> &seq, types::str const &);

    DEFINE_FUNCTOR(pythonic::array::array, frombytes);
  } // namespace array
} // namespace array
PYTHONIC_NS_END
#endif
