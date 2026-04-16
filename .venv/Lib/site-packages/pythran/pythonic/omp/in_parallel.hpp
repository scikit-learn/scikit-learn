#ifndef PYTHONIC_OMP_IN_PARALLEL_HPP
#define PYTHONIC_OMP_IN_PARALLEL_HPP

#include "pythonic/include/omp/in_parallel.hpp"

#include "pythonic/utils/functor.hpp"
#include <omp.h>

PYTHONIC_NS_BEGIN

namespace omp
{

  bool in_parallel()
  {
    return omp_in_parallel();
  }
} // namespace omp
PYTHONIC_NS_END

#endif
