#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_READLINES_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_READLINES_HPP

#include "pythonic/include/builtins/file/readlines.hpp"

PYTHONIC_NS_BEGIN
namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(readlines, builtins::file::functor::readlines);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
