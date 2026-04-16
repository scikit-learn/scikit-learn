#ifndef UFUNC_INAME
#error missing UFUNC_INAME
#endif

// clang-format off
#include INCLUDE_FILE(pythonic/operator_,UFUNC_INAME)
// clang-format on
#include "pythonic/numpy/reduce.hpp"
#include "pythonic/utils/functor.hpp"
