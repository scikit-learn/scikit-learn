#ifndef PYTHONIC_INCLUDE_OMP_SET_NESTED_HPP
#define PYTHONIC_INCLUDE_OMP_SET_NESTED_HPP

#include "pythonic/include/utils/functor.hpp"
#include <omp.h>

PYTHONIC_NS_BEGIN

namespace omp
{

  void set_nested(long val);

  DEFINE_FUNCTOR(pythonic::omp, set_nested);
} // namespace omp
PYTHONIC_NS_END

#endif
