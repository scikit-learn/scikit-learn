#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_FLUSH_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_FLUSH_HPP

#include "pythonic/include/builtins/file/flush.hpp"

PYTHONIC_NS_BEGIN

namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(flush, builtins::file::functor::flush);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
