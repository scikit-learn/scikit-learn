#ifndef PYTHONIC_NUMPY_BINARYREPR_HPP
#define PYTHONIC_NUMPY_BINARYREPR_HPP

#include "pythonic/include/numpy/binary_repr.hpp"

#include "pythonic/numpy/base_repr.hpp"
#include "pythonic/utils/allocate.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline char *int2bin(long a, char *buffer, int buf_size)
    {
      buffer += (buf_size - 1);
      for (int i = 0; i < buf_size; ++i) {
        *buffer-- = (a & 1) + '0';
        a >>= 1;
      }
      return buffer;
    }
  } // namespace details

  inline types::str binary_repr(long number, types::none_type width)
  {
    return base_repr(number, 2);
  }

  inline types::str binary_repr(long number, long width)
  {
    types::str out = binary_repr(std::abs(number));
    if (number >= 0)
      return base_repr(number, 2, width - out.size());
    else {
      char *mem = utils::allocate<char>(width);
      details::int2bin(number, mem, width);
      auto res = types::str(mem, width);
      return res;
    }
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
