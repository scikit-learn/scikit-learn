#ifndef PYTHONIC_UTILS_ALLOCATE_HPP
#define PYTHONIC_UTILS_ALLOCATE_HPP

#include "pythonic/include/utils/allocate.hpp"
PYTHONIC_NS_BEGIN

namespace utils
{
#ifdef PYTHRAN_TRACE_ALLOCATION
  size_t pythran_allocation_site;
#endif
} // namespace utils

PYTHONIC_NS_END

#endif
