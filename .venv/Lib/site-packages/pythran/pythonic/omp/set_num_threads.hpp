#ifndef PYTHONIC_OMP_SET_NUM_THREADS_HPP
#define PYTHONIC_OMP_SET_NUM_THREADS_HPP

#include "pythonic/include/omp/set_num_threads.hpp"

#include "pythonic/utils/functor.hpp"
#include <omp.h>

PYTHONIC_NS_BEGIN

namespace omp
{
  void set_num_threads(long num_threads)
  {
    return omp_set_num_threads(num_threads);
  }
} // namespace omp
PYTHONIC_NS_END

#endif
