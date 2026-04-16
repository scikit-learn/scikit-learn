#ifndef PYTHONIC_OMP_GET_WTIME_HPP
#define PYTHONIC_OMP_GET_WTIME_HPP

#include "pythonic/include/omp/get_wtime.hpp"

#include "pythonic/utils/functor.hpp"
#include <omp.h>

PYTHONIC_NS_BEGIN

namespace omp
{

  long get_wtime()
  {
    return omp_get_wtime();
  }
} // namespace omp
PYTHONIC_NS_END

#endif
