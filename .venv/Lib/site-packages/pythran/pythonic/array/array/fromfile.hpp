#ifndef PYTHONIC_ARRAY_ARRAY_FROMFILE_HPP
#define PYTHONIC_ARRAY_ARRAY_FROMFILE_HPP

#include "pythonic/include/types/array.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace array
{

  namespace array
  {

    template <class T>
    types::none_type fromfile(types::array<T> &seq, types::file &f, long n)
    {
      long p = seq.size();
      seq.resize(p + n);
      f.read_as(n, seq.data() + p);
      return {};
    }

  } // namespace array
} // namespace array
PYTHONIC_NS_END
#endif
