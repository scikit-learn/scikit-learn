#ifndef PYTHONIC_NUMPY_INT64_HPP
#define PYTHONIC_NUMPY_INT64_HPP

#include "pythonic/include/numpy/int64.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline int64_t int64()
    {
      return int64_t();
    }

    template <class V>
    int64_t int64(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME int64
#define NUMPY_NARY_FUNC_SYM details::int64
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::int64>::convert(numpy::functor::int64 const &c)
{
  return (PyObject *)&PyInt64ArrType_Type;
}

inline bool from_python<numpy::functor::int64>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyInt64ArrType_Type;
}

inline numpy::functor::int64 from_python<numpy::functor::int64>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
