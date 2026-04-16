#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_WRITELINES_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_WRITELINES_HPP

#include "pythonic/include/builtins/file/writelines.hpp"

PYTHONIC_NS_BEGIN
namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(writelines, builtins::file::functor::writelines);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
