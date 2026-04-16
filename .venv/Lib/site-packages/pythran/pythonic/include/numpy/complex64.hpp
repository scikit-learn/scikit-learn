#ifndef PYTHONIC_INCLUDE_NUMPY_COMPLEX64_HPP
#define PYTHONIC_INCLUDE_NUMPY_COMPLEX64_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/types/numpy_op_helper.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/meta.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {
    std::complex<float> complex64();
    template <class V>
    std::complex<float> complex64(V v);
    template <class T>
    std::complex<float> complex64(std::complex<T> v);
  } // namespace details

#define NUMPY_NARY_FUNC_NAME complex64
#define NUMPY_NARY_FUNC_SYM details::complex64
#define NUMPY_NARY_EXTRA_METHOD using type = std::complex<float>;
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<numpy::functor::complex64> {
  static PyObject *convert(numpy::functor::complex64 const &c);
};

template <>
struct from_python<numpy::functor::complex64> {
  static bool is_convertible(PyObject *obj);
  static numpy::functor::complex64 convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
