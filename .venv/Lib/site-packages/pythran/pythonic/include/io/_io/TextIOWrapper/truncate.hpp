#ifndef PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_TRUNCATE_HPP
#define PYTHONIC_INCLUDE_IO__IO_TEXTIOWRAPPER_TRUNCATE_HPP

#include "pythonic/include/builtins/file/truncate.hpp"

PYTHONIC_NS_BEGIN
namespace io
{

  namespace _io
  {
    namespace TextIOWrapper
    {
      USING_FUNCTOR(truncate, builtins::file::functor::truncate);
    }
  } // namespace _io
} // namespace io
PYTHONIC_NS_END
#endif
