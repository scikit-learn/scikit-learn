#ifndef PYTHONIC_ARRAY_ARRAY_FROMBYTES_HPP
#define PYTHONIC_ARRAY_ARRAY_FROMBYTES_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/types/array.hpp"
#include "pythonic/types/str.hpp"

PYTHONIC_NS_BEGIN

namespace array
{

  namespace array
  {

    template <class T>
    types::none_type frombytes(types::array<T> &seq, types::str const &s)
    {
      long size = seq.size();
      seq.resize(size + s.size() / sizeof(T));
      memcpy(seq.data() + size, s.c_str(), s.size());
      return {};
    }

  } // namespace array
} // namespace array
PYTHONIC_NS_END
#endif
