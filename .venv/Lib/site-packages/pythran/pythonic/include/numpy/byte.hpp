#ifndef PYTHONIC_INCLUDE_NUMPY_BYTE_HPP
#define PYTHONIC_INCLUDE_NUMPY_BYTE_HPP

#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    char byte();
    template <class V>
    char byte(V v);
  } // namespace details

#define NUMPY_NARY_FUNC_NAME byte
#define NUMPY_NARY_FUNC_SYM details::byte
#define NUMPY_NARY_EXTRA_METHOD using type = char;
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<numpy::functor::byte> {
  static PyObject *convert(numpy::functor::byte const &c);
};

template <>
struct from_python<numpy::functor::byte> {
  static bool is_convertible(PyObject *obj);
  static numpy::functor::byte convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
