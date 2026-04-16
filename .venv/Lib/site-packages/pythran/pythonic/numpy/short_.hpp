#ifndef PYTHONIC_NUMPY_SHORT__HPP
#define PYTHONIC_NUMPY_SHORT__HPP

#include "pythonic/include/numpy/short_.hpp"

#include "pythonic/types/numpy_op_helper.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/meta.hpp"
#include "pythonic/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    inline short short_()
    {
      return {};
    }

    template <class V>
    short short_(V v)
    {
      return v;
    }
  } // namespace details

#define NUMPY_NARY_FUNC_NAME short_
#define NUMPY_NARY_FUNC_SYM details::short_
#include "pythonic/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayscalars.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<numpy::functor::short_>::convert(numpy::functor::short_ const &c)
{
  return (PyObject *)&PyShortArrType_Type;
}

inline bool from_python<numpy::functor::short_>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyShortArrType_Type;
}

inline numpy::functor::short_ from_python<numpy::functor::short_>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
