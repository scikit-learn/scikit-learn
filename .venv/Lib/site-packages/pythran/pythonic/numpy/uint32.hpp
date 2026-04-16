#ifndef PYTHONIC_NUMPY_UINT32_HPP
#define PYTHONIC_NUMPY_UINT32_HPP

#include "pythonic/include/numpy/uint32.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline uint32_t uint32()
    {
      return uint32_t();
    }

    template <class V>
    uint32_t uint32(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME uint32
#define NUMPY_NARY_FUNC_SYM details::uint32
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::uint32>::convert(numpy::functor::uint32 const &c)
{
  return (PyObject *)&PyUInt32ArrType_Type;
}

inline bool from_python<numpy::functor::uint32>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyUInt32ArrType_Type;
}

inline numpy::functor::uint32 from_python<numpy::functor::uint32>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
