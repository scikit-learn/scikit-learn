#ifndef PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_KWONLY_HPP
#define PYTHONIC_INCLUDE_BUILTIN_PYTHRAN_KWONLY_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace types
{
  struct kwonly {
  };
} // namespace types

namespace builtins
{
  namespace pythran
  {
    inline types::kwonly kwonly()
    {
      return {};
    };

    DEFINE_FUNCTOR(pythonic::builtins::pythran, kwonly);
  } // namespace pythran
} // namespace builtins
PYTHONIC_NS_END

#endif
