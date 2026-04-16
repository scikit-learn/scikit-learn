#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_TELL_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_TELL_HPP

#include "pythonic/include/builtins/file/tell.hpp"

PYTHONIC_NS_BEGIN
namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(tell, builtins::file::functor::tell);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
