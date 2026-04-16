#ifndef PYTHONIC_NUMPY_FLOAT128_HPP
#define PYTHONIC_NUMPY_FLOAT128_HPP

#include "pythonic/include/numpy/float128.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {

    template <class V>
    long double float128(V v)
    {
      return static_cast<long double>(v);
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME float128
#define NUMPY_NARY_FUNC_SYM details::float128
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::float128>::convert(numpy::functor::float128 const &c)
{
  return (PyObject *)&PyLongDoubleArrType_Type;
}

inline bool from_python<numpy::functor::float128>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyLongDoubleArrType_Type;
}

inline numpy::functor::float128 from_python<numpy::functor::float128>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
