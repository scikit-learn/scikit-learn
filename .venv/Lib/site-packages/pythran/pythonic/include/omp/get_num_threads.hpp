#ifndef PYTHONIC_INCLUDE_OMP_GET_NUM_THREADS_HPP
#define PYTHONIC_INCLUDE_OMP_GET_NUM_THREADS_HPP

#include <omp.h>

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace omp
{

  long get_num_threads();

  DEFINE_FUNCTOR(pythonic::omp, get_num_threads);
} // namespace omp
PYTHONIC_NS_END

#endif
