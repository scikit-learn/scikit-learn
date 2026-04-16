#ifndef PYTHONIC_NUMPY_COMPLEX128_HPP
#define PYTHONIC_NUMPY_COMPLEX128_HPP

#include "pythonic/include/numpy/complex128.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {

    inline std::complex<double> complex128()
    {
      return {};
    }

    template <class V>
    std::complex<double> complex128(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME complex128
#define NUMPY_NARY_FUNC_SYM details::complex128
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::complex128>::convert(numpy::functor::complex128 const &c)
{
  return (PyObject *)&PyCDoubleArrType_Type;
}

inline bool from_python<numpy::functor::complex128>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyCDoubleArrType_Type;
}

inline numpy::functor::complex128 from_python<numpy::functor::complex128>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
