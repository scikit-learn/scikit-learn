#ifndef PYTHONIC_NUMPY_COMPLEX64_HPP
#define PYTHONIC_NUMPY_COMPLEX64_HPP

#include "pythonic/include/numpy/complex64.hpp"

#include "pythonic/types/complex.hpp"
#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {

    inline std::complex<float> complex64()
    {
      return {};
    }

    template <class V>
    std::complex<float> complex64(V v)
    {
      return v;
    }

    template <class T>
    std::complex<float> complex64(std::complex<T> v)
    {
      return {(float)v.real(), (float)v.imag()};
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME complex64
#define NUMPY_NARY_FUNC_SYM details::complex64
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::complex64>::convert(numpy::functor::complex64 const &c)
{
  return (PyObject *)&PyCFloatArrType_Type;
}

inline bool from_python<numpy::functor::complex64>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyCFloatArrType_Type;
}

inline numpy::functor::complex64 from_python<numpy::functor::complex64>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
