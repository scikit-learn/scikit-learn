
#ifndef PYTHONIC_CORE_HPP
#define PYTHONIC_CORE_HPP

#define PYTHONIC_NS_BEGIN                                                                          \
  namespace                                                                                        \
  {                                                                                                \
    namespace pythonic                                                                             \
    {
#define PYTHONIC_NS_END                                                                            \
  }                                                                                                \
  }

// mostly to flag '_' as unused in generated code
#ifdef WIN32
#define PYTHRAN_UNUSED
#else
#define PYTHRAN_UNUSED __attribute__((unused))
#endif

// for backward compatibility
#ifdef USE_BOOST_SIMD
#define USE_XSIMD
#endif

#define STR_(M) #M
// clang-format off
#define INCLUDE_FILE(U, M) STR_(U/M.hpp)
// clang-format on

#ifdef ENABLE_PYTHON_MODULE

#if defined(Py_LIMITED_API) && Py_LIMITED_API < 0x03070000
#error Pythran only supports Py_LIMITED_API from 0x03070000 upward
#endif

// Define python's visibility macros
#include "pyconfig.h"

// Some version of python define that macro on Windows, and it breaks
// compilation of some C++ headers.
#ifdef copysign
#undef copysign
#endif
#endif

#include "pythonic/types/assignable.hpp"
#include "pythonic/types/combined.hpp"

#include "pythonic/types/float.hpp"
#include "pythonic/types/int.hpp"
#include "pythonic/types/slice.hpp"

#endif
