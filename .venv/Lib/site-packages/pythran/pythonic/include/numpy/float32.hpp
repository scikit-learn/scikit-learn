#ifndef PYTHONIC_INCLUDE_NUMPY_FLOAT32_HPP
#define PYTHONIC_INCLUDE_NUMPY_FLOAT32_HPP

#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    float float32();
    template <class V>
    float float32(V v);
  } // namespace details

#define NUMPY_NARY_FUNC_NAME float32
#define NUMPY_NARY_FUNC_SYM details::float32
#define NUMPY_NARY_EXTRA_METHOD using type = float;
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<numpy::functor::float32> {
  static PyObject *convert(numpy::functor::float32 const &c);
};

template <>
struct from_python<numpy::functor::float32> {
  static bool is_convertible(PyObject *obj);
  static numpy::functor::float32 convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
