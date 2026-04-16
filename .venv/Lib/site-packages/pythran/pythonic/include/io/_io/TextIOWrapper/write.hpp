#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_WRITE_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_WRITE_HPP

#include "pythonic/include/builtins/file/write.hpp"

PYTHONIC_NS_BEGIN
namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(write, builtins::file::functor::write);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
