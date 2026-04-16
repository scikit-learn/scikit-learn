#ifndef PYTHONIC_NUMPY_UINT16_HPP
#define PYTHONIC_NUMPY_UINT16_HPP

#include "pythonic/include/numpy/uint16.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline uint16_t uint16()
    {
      return uint16_t();
    }

    template <class V>
    uint16_t uint16(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME uint16
#define NUMPY_NARY_FUNC_SYM details::uint16
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::uint16>::convert(numpy::functor::uint16 const &c)
{
  return (PyObject *)&PyUInt16ArrType_Type;
}

inline bool from_python<numpy::functor::uint16>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyUInt16ArrType_Type;
}

inline numpy::functor::uint16 from_python<numpy::functor::uint16>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
