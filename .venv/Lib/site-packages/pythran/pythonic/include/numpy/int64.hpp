#ifndef PYTHONIC_INCLUDE_NUMPY_INT64_HPP
#define PYTHONIC_INCLUDE_NUMPY_INT64_HPP

#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    int64_t int64();
    template <class V>
    int64_t int64(V v);
  } // namespace details

#define NUMPY_NARY_FUNC_NAME int64
#define NUMPY_NARY_FUNC_SYM details::int64
#define NUMPY_NARY_EXTRA_METHOD using type = int64_t;
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<numpy::functor::int64> {
  static PyObject *convert(numpy::functor::int64 const &c);
};

template <>
struct from_python<numpy::functor::int64> {
  static bool is_convertible(PyObject *obj);
  static numpy::functor::int64 convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
