#ifndef PYTHONIC_NUMPY_UINT8_HPP
#define PYTHONIC_NUMPY_UINT8_HPP

#include "pythonic/include/numpy/uint8.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline uint8_t uint8()
    {
      return uint8_t();
    }

    template <class V>
    uint8_t uint8(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME uint8
#define NUMPY_NARY_FUNC_SYM details::uint8
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::uint8>::convert(numpy::functor::uint8 const &c)
{
  return (PyObject *)&PyUInt8ArrType_Type;
}

inline bool from_python<numpy::functor::uint8>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyUInt8ArrType_Type;
}

inline numpy::functor::uint8 from_python<numpy::functor::uint8>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
