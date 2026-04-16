#ifndef PYTHONIC_INCLUDE_NUMPY_UBYTE_HPP
#define PYTHONIC_INCLUDE_NUMPY_UBYTE_HPP

#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    unsigned char ubyte();
    template <class V>
    unsigned char ubyte(V v);
  } // namespace details

#define NUMPY_NARY_FUNC_NAME ubyte
#define NUMPY_NARY_FUNC_SYM details::ubyte
#define NUMPY_NARY_EXTRA_METHOD using type = unsigned char;
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<numpy::functor::ubyte> {
  static PyObject *convert(numpy::functor::ubyte const &c);
};

template <>
struct from_python<numpy::functor::ubyte> {
  static bool is_convertible(PyObject *obj);
  static numpy::functor::ubyte convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
