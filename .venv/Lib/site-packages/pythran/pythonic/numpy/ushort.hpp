#ifndef PYTHONIC_NUMPY_USHORT_HPP
#define PYTHONIC_NUMPY_USHORT_HPP

#include "pythonic/include/numpy/ushort.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline unsigned short ushort()
    {
      return {};
    }

    template <class V>
    unsigned short ushort(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME ushort
#define NUMPY_NARY_FUNC_SYM details::ushort
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::ushort>::convert(numpy::functor::ushort const &c)
{
  return (PyObject *)&PyUShortArrType_Type;
}

inline bool from_python<numpy::functor::ushort>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyUShortArrType_Type;
}

inline numpy::functor::ushort from_python<numpy::functor::ushort>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
