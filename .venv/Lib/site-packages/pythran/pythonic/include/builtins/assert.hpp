#ifndef PYTHONIC_INCLUDE_BUILTIN_ASSERT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_ASSERT_HPP

#include "pythonic/include/types/str.hpp"

PYTHONIC_NS_BEGIN
void pythran_assert(bool cond);
void pythran_assert(bool cond, types::str const &what);
PYTHONIC_NS_END

#endif
