#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_NEXT_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_NEXT_HPP

#include "pythonic/include/builtins/file/next.hpp"

PYTHONIC_NS_BEGIN
namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(next, builtins::file::functor::next);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
