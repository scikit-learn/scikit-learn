#ifndef PYTHONIC_NUMPY_UINTC_HPP
#define PYTHONIC_NUMPY_UINTC_HPP

#include "pythonic/include/numpy/uintc.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline unsigned uintc()
    {
      return {};
    }

    template <class V>
    unsigned uintc(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME uintc
#define NUMPY_NARY_FUNC_SYM details::uintc
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::uintc>::convert(numpy::functor::uintc const &c)
{
  if (sizeof(unsigned) == 4)
    return (PyObject *)&PyUInt32ArrType_Type;
  else
    return (PyObject *)&PyUInt64ArrType_Type;
}

inline bool from_python<numpy::functor::uintc>::is_convertible(PyObject *obj)
{
  if (sizeof(unsigned) == 4)
    return obj == (PyObject *)&PyUInt32ArrType_Type;
  else
    return obj == (PyObject *)&PyUInt64ArrType_Type;
}

inline numpy::functor::uintc from_python<numpy::functor::uintc>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
