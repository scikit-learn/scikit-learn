#ifndef PYTHONIC_INCLUDE_NUMPY_FLOAT64_HPP
#define PYTHONIC_INCLUDE_NUMPY_FLOAT64_HPP

#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {

    double float64();
    template <class V>
    double float64(V v);
  } // namespace details

#define NUMPY_NARY_FUNC_NAME float64
#define NUMPY_NARY_FUNC_SYM details::float64
#define NUMPY_NARY_EXTRA_METHOD using type = double;
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<numpy::functor::float64> {
  static PyObject *convert(numpy::functor::float64 const &c);
};

template <>
struct from_python<numpy::functor::float64> {
  static bool is_convertible(PyObject *obj);
  static numpy::functor::float64 convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
