#ifndef PYTHONIC_INCLUDE_NUMPY_UINT8_HPP
#define PYTHONIC_INCLUDE_NUMPY_UINT8_HPP

#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace details
  {

    uint8_t uint8();
    template <class V>
    uint8_t uint8(V v);
  } // namespace details

#define NUMPY_NARY_FUNC_NAME uint8
#define NUMPY_NARY_FUNC_SYM details::uint8
#define NUMPY_NARY_EXTRA_METHOD using type = uint8_t;
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<numpy::functor::uint8> {
  static PyObject *convert(numpy::functor::uint8 const &c);
};

template <>
struct from_python<numpy::functor::uint8> {
  static bool is_convertible(PyObject *obj);
  static numpy::functor::uint8 convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
