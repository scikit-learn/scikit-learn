#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_READ_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_READ_HPP

#include "pythonic/include/builtins/file/read.hpp"

PYTHONIC_NS_BEGIN
namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(read, builtins::file::functor::read);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
