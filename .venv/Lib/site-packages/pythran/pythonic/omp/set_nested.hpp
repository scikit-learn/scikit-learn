#ifndef PYTHONIC_OMP_SET_NESTED_HPP
#define PYTHONIC_OMP_SET_NESTED_HPP

#include "pythonic/include/omp/set_nested.hpp"

#include "pythonic/utils/functor.hpp"
#include <omp.h>

PYTHONIC_NS_BEGIN

namespace omp
{

  void set_nested(long val)
  {
    return omp_set_nested(val);
  }
} // namespace omp
PYTHONIC_NS_END

#endif
