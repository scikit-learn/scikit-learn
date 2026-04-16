#ifndef PYTHONIC_NUMPY_INTC_HPP
#define PYTHONIC_NUMPY_INTC_HPP

#include "pythonic/include/numpy/intc.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline int intc()
    {
      return {};
    }

    template <class V>
    int intc(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME intc
#define NUMPY_NARY_FUNC_SYM details::intc
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::intc>::convert(numpy::functor::intc const &c)
{
  if (sizeof(int) == 4)
    return (PyObject *)&PyInt32ArrType_Type;
  else
    return (PyObject *)&PyInt64ArrType_Type;
}

inline bool from_python<numpy::functor::intc>::is_convertible(PyObject *obj)
{
  if (sizeof(int) == 4)
    return obj == (PyObject *)&PyInt32ArrType_Type;
  else
    return obj == (PyObject *)&PyInt64ArrType_Type;
}

inline numpy::functor::intc from_python<numpy::functor::intc>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
