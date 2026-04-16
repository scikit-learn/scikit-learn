#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_CLOSE_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_CLOSE_HPP

#include "pythonic/include/builtins/file/close.hpp"

PYTHONIC_NS_BEGIN

namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(close, builtins::file::functor::close);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
