#ifndef PYTHONIC_NUMPY_INT32_HPP
#define PYTHONIC_NUMPY_INT32_HPP

#include "pythonic/include/numpy/int32.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline int32_t int32()
    {
      return int32_t();
    }

    template <class V>
    int32_t int32(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME int32
#define NUMPY_NARY_FUNC_SYM details::int32
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::int32>::convert(numpy::functor::int32 const &c)
{
  return (PyObject *)&PyInt32ArrType_Type;
}

inline bool from_python<numpy::functor::int32>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyInt32ArrType_Type;
}

inline numpy::functor::int32 from_python<numpy::functor::int32>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
