#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_READLINE_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_READLINE_HPP

#include "pythonic/include/builtins/file/readline.hpp"

PYTHONIC_NS_BEGIN
namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(readline, builtins::file::functor::readline);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
