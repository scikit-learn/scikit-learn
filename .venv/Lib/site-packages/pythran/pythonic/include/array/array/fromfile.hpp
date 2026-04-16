#ifndef PYTHONIC_INCLUDE_ARRAY_ARRAY_FROMFILE_HPP
#define PYTHONIC_INCLUDE_ARRAY_ARRAY_FROMFILE_HPP

#include "pythonic/include/types/array.hpp"
#include "pythonic/include/types/file.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace array
{

  namespace array
  {

    template <class T>
    types::none_type fromfile(types::array<T> &seq, types::file &f, long n);

    template <class T>
    types::none_type fromfile(types::array<T> &seq, types::file &&f, long n)
    {
      return fromfile(seq, f, n);
    }

    template <class T>
    types::none_type fromfile(types::array<T> &&seq, types::file &&f, long n)
    {
      return fromfile(seq, f, n);
    }

    template <class T>
    types::none_type fromfile(types::array<T> &&seq, types::file &f, long n)
    {
      return fromfile(seq, f, n);
    }

    DEFINE_FUNCTOR(pythonic::array::array, fromfile);
  } // namespace array
} // namespace array
PYTHONIC_NS_END
#endif
