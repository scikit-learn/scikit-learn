#ifndef PYTHONIC_OMP_GET_NUM_THREADS_HPP
#define PYTHONIC_OMP_GET_NUM_THREADS_HPP

#include "pythonic/include/omp/get_num_threads.hpp"

#include "pythonic/utils/functor.hpp"
#include <omp.h>

PYTHONIC_NS_BEGIN

namespace omp
{
  long get_num_threads()
  {
    return omp_get_num_threads();
  }
} // namespace omp
PYTHONIC_NS_END

#endif
