#ifndef PYTHONIC_NUMPY_UINT_HPP
#define PYTHONIC_NUMPY_UINT_HPP

#include "pythonic/include/numpy/uint.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline unsigned long uint()
    {
      return {};
    }

    template <class V>
    unsigned long uint(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME uint
#define NUMPY_NARY_FUNC_SYM details::uint
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::uint>::convert(numpy::functor::uint const &c)
{
  return (PyObject *)&PyUIntArrType_Type;
}

inline bool from_python<numpy::functor::uint>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyUIntArrType_Type ||
#if NPY_SIZEOF_INTP == NPY_SIZEOF_LONG
         obj == (PyObject *)&PyUInt64ArrType_Type
#else
         obj == (PyObject *)&PyUInt32ArrType_Type
#endif
      ;
}

inline numpy::functor::uint from_python<numpy::functor::uint>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
