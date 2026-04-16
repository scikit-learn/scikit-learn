#ifndef PYTHONIC_INCLUDE_OMP_SET_NUM_THREADS_HPP
#define PYTHONIC_INCLUDE_OMP_SET_NUM_THREADS_HPP

#include <omp.h>

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace omp
{

  void set_num_threads(long);

  DEFINE_FUNCTOR(pythonic::omp, set_num_threads);
} // namespace omp
PYTHONIC_NS_END

#endif
