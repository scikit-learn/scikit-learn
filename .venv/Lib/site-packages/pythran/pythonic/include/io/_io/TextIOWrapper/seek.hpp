#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_SEEK_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_SEEK_HPP

#include "pythonic/include/builtins/file/seek.hpp"

PYTHONIC_NS_BEGIN
namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(seek, builtins::file::functor::seek);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
