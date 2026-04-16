#ifndef PYTHONIC_BUILTIN_COMPLEX_HPP
#define PYTHONIC_BUILTIN_COMPLEX_HPP

#include "pythonic/include/builtins/complex.hpp"

#include "pythonic/types/complex.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace functor
  {
    inline complex::type complex::operator()(double v0, double v1) const
    {
      return {v0, v1};
    }
  } // namespace functor
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<builtins::functor::complex>::convert(builtins::functor::complex const &c)
{
  return (PyObject *)&PyComplex_Type;
}

inline bool from_python<builtins::functor::complex>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyComplex_Type || obj == (PyObject *)&PyCDoubleArrType_Type;
}

inline builtins::functor::complex from_python<builtins::functor::complex>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
