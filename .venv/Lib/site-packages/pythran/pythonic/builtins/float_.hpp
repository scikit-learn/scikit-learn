#ifndef PYTHONIC_BUILTIN_FLOAT_HPP
#define PYTHONIC_BUILTIN_FLOAT_HPP

#include "pythonic/include/builtins/float_.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace functor
  {
    template <class T>
    float_::type float_::operator()(T &&t) const
    {
      return static_cast<float_::type>(t);
    }

    inline float_::type float_::operator()() const
    {
      return 0.;
    }
  } // namespace functor
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<builtins::functor::float_>::convert(builtins::functor::float_ const &c)
{
  return (PyObject *)&PyFloat_Type;
}

inline bool from_python<builtins::functor::float_>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyFloat_Type || obj == (PyObject *)&PyDoubleArrType_Type;
}

inline builtins::functor::float_ from_python<builtins::functor::float_>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
