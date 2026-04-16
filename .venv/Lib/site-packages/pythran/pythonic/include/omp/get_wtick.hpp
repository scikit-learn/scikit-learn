#ifndef PYTHONIC_INCLUDE_OMP_GET_WTICK_HPP
#define PYTHONIC_INCLUDE_OMP_GET_WTICK_HPP

#include "pythonic/include/utils/functor.hpp"
#include <omp.h>

PYTHONIC_NS_BEGIN

namespace omp
{
  long get_wtick();

  DEFINE_FUNCTOR(pythonic::omp, get_wtick);
} // namespace omp
PYTHONIC_NS_END

#endif
