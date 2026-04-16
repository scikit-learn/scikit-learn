#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_ISATTY_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_ISATTY_HPP

#include "pythonic/include/builtins/file/isatty.hpp"

PYTHONIC_NS_BEGIN

namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(isatty, builtins::file::functor::isatty);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
