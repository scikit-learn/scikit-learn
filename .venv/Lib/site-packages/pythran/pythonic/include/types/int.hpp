#ifndef PYTHONIC_INCLUDE_TYPES_INT_HPP
#define PYTHONIC_INCLUDE_TYPES_INT_HPP

#include "pythonic/include/types/attr.hpp"

PYTHONIC_NS_BEGIN
namespace builtins
{
  template <class T>
  typename std::enable_if<std::is_integral<T>::value, T>::value getattr(types::attr::REAL, T self);
  template <class T>
  typename std::enable_if<std::is_integral<T>::value, T>::value getattr(types::attr::IMAG, T self);
} // namespace builtins
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

#define PYTHONIC_INT_TO_PYTHON(TYPE)                                                               \
  template <>                                                                                      \
  struct to_python<TYPE> {                                                                         \
    static PyObject *convert(TYPE l);                                                              \
  }

PYTHONIC_INT_TO_PYTHON(char);
PYTHONIC_INT_TO_PYTHON(unsigned char);
PYTHONIC_INT_TO_PYTHON(signed char);
PYTHONIC_INT_TO_PYTHON(unsigned short);
PYTHONIC_INT_TO_PYTHON(signed short);
PYTHONIC_INT_TO_PYTHON(unsigned int);
PYTHONIC_INT_TO_PYTHON(signed int);
PYTHONIC_INT_TO_PYTHON(unsigned long);
PYTHONIC_INT_TO_PYTHON(signed long);
PYTHONIC_INT_TO_PYTHON(unsigned long long);
PYTHONIC_INT_TO_PYTHON(signed long long);

#undef PYTHONIC_INT_TO_PYTHON

#define PYTHONIC_INT_FROM_PYTHON(TYPE)                                                             \
  template <>                                                                                      \
  struct from_python<TYPE> {                                                                       \
    static bool is_convertible(PyObject *obj);                                                     \
    static TYPE convert(PyObject *obj);                                                            \
  }

PYTHONIC_INT_FROM_PYTHON(unsigned char);
PYTHONIC_INT_FROM_PYTHON(signed char);
PYTHONIC_INT_FROM_PYTHON(unsigned short);
PYTHONIC_INT_FROM_PYTHON(signed short);
PYTHONIC_INT_FROM_PYTHON(unsigned int);
PYTHONIC_INT_FROM_PYTHON(signed int);
PYTHONIC_INT_FROM_PYTHON(unsigned long);
PYTHONIC_INT_FROM_PYTHON(signed long);
PYTHONIC_INT_FROM_PYTHON(unsigned long long);
PYTHONIC_INT_FROM_PYTHON(signed long long);

#undef PYTHONIC_INT_FROM_PYTHON

PYTHONIC_NS_END
#endif

#endif
