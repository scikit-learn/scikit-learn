#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_FILENO_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_FILENO_HPP

#include "pythonic/include/builtins/file/fileno.hpp"

PYTHONIC_NS_BEGIN

namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(fileno, builtins::file::functor::fileno);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
