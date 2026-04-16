#ifndef PYTHONIC_BUILTIN_ASSERT_HPP
#define PYTHONIC_BUILTIN_ASSERT_HPP

#include "pythonic/include/builtins/assert.hpp"

#include "pythonic/builtins/AssertionError.hpp"
#include "pythonic/types/str.hpp"

PYTHONIC_NS_BEGIN

inline void pythran_assert(bool cond)
{
#ifndef NDEBUG
  if (!cond)
    throw types::AssertionError();
#endif
}

inline void pythran_assert(bool cond, types::str const &what)
{
#ifndef NDEBUG
  if (!cond)
    throw types::AssertionError(what);
#endif
}
PYTHONIC_NS_END

#endif
