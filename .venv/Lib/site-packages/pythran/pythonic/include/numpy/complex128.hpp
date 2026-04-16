#ifndef PYTHONIC_INCLUDE_NUMPY_COMPLEX128_HPP
#define PYTHONIC_INCLUDE_NUMPY_COMPLEX128_HPP

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
    std::complex<double> complex128();
    template <class V>
    std::complex<double> complex128(V v);
  } // namespace details

#define NUMPY_NARY_FUNC_NAME complex128
#define NUMPY_NARY_FUNC_SYM details::complex128
#define NUMPY_NARY_EXTRA_METHOD using type = std::complex<double>;
#include "pythonic/include/types/numpy_nary_expr.hpp"
} // namespace numpy
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<numpy::functor::complex128> {
  static PyObject *convert(numpy::functor::complex128 const &c);
};

template <>
struct from_python<numpy::functor::complex128> {
  static bool is_convertible(PyObject *obj);
  static numpy::functor::complex128 convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
