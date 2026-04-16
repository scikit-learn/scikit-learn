#ifndef PYTHONIC_INCLUDE_BUILTIN_COMPLEX_HPP
#define PYTHONIC_INCLUDE_BUILTIN_COMPLEX_HPP

#include "pythonic/include/types/complex.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace functor
  {
    struct complex {
      using callable = void;
      using type = std::complex<double>;
      // TODO: doesn't handle string as first argument
      type operator()(double v0 = 0, double v1 = 0) const;
      friend std::ostream &operator<<(std::ostream &os, complex)
      {
        return os << "complex";
      }
    };
  } // namespace functor
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<builtins::functor::complex> {
  static PyObject *convert(builtins::functor::complex const &c);
};

template <>
struct from_python<builtins::functor::complex> {
  static bool is_convertible(PyObject *obj);
  static builtins::functor::complex convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
