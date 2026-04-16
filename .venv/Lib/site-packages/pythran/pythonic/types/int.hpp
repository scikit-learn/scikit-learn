#ifndef PYTHONIC_TYPES_INT_HPP
#define PYTHONIC_TYPES_INT_HPP
#include <iostream>

#include "pythonic/include/types/int.hpp"
#include "pythonic/types/attr.hpp"

PYTHONIC_NS_BEGIN
namespace builtins
{
  template <class T>
  std::enable_if_t<std::is_integral<T>::value, T> getattr(types::attr::REAL, T self)
  {
    return self;
  }
  template <class T>
  std::enable_if_t<std::is_integral<T>::value, T> getattr(types::attr::IMAG, T self)
  {
    return T(0);
  }
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "numpy/arrayobject.h"
#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <class T>
struct c_type_to_numpy_type : c_type_to_numpy_type<decltype(std::declval<T>()())> {
};

template <>
struct c_type_to_numpy_type<long double> : std::integral_constant<int, NPY_LONGDOUBLE> {
};

template <>
struct c_type_to_numpy_type<double> : std::integral_constant<int, NPY_DOUBLE> {
};

template <>
struct c_type_to_numpy_type<float> : std::integral_constant<int, NPY_FLOAT> {
};

template <>
struct c_type_to_numpy_type<std::complex<float>> : std::integral_constant<int, NPY_CFLOAT> {
};

template <>
struct c_type_to_numpy_type<std::complex<double>> : std::integral_constant<int, NPY_CDOUBLE> {
};

template <>
struct c_type_to_numpy_type<std::complex<long double>>
    : std::integral_constant<int, NPY_CLONGDOUBLE> {
};

template <>
struct c_type_to_numpy_type<signed long long> : std::integral_constant<int, NPY_LONGLONG> {
};

template <>
struct c_type_to_numpy_type<unsigned long long> : std::integral_constant<int, NPY_ULONGLONG> {
};

template <>
struct c_type_to_numpy_type<signed long> : std::integral_constant<int, NPY_LONG> {
};

template <>
struct c_type_to_numpy_type<unsigned long> : std::integral_constant<int, NPY_ULONG> {
};

template <>
struct c_type_to_numpy_type<signed int> : std::integral_constant<int, NPY_INT> {
};

template <>
struct c_type_to_numpy_type<unsigned int> : std::integral_constant<int, NPY_UINT> {
};

template <>
struct c_type_to_numpy_type<signed short> : std::integral_constant<int, NPY_SHORT> {
};

template <>
struct c_type_to_numpy_type<unsigned short> : std::integral_constant<int, NPY_USHORT> {
};

template <>
struct c_type_to_numpy_type<char> : std::integral_constant<int, NPY_BYTE> {
};

template <>
struct c_type_to_numpy_type<signed char> : std::integral_constant<int, NPY_BYTE> {
};

template <>
struct c_type_to_numpy_type<unsigned char> : std::integral_constant<int, NPY_UBYTE> {
};

template <>
struct c_type_to_numpy_type<bool> : std::integral_constant<int, NPY_BOOL> {
};

#ifndef PyInt_FromLong
#define PyInt_FromLong PyLong_FromLong

#ifndef PyInt_CheckExact
#define PyInt_CheckExact PyLong_CheckExact
#endif
#ifndef PyInt_AsLong
#define PyInt_AsLong PyLong_AsLong
#endif
#endif

#define PYTHONIC_INT_TO_PYTHON(TYPE)                                                               \
  inline PyObject *to_python<TYPE>::convert(TYPE l)                                                \
  {                                                                                                \
    return PyArray_Scalar(&l, PyArray_DescrFromType(c_type_to_numpy_type<TYPE>::value), nullptr);  \
  }

PYTHONIC_INT_TO_PYTHON(char)
PYTHONIC_INT_TO_PYTHON(unsigned char)
PYTHONIC_INT_TO_PYTHON(signed char)
PYTHONIC_INT_TO_PYTHON(unsigned short)
PYTHONIC_INT_TO_PYTHON(signed short)
PYTHONIC_INT_TO_PYTHON(unsigned int)
PYTHONIC_INT_TO_PYTHON(signed int)
PYTHONIC_INT_TO_PYTHON(unsigned long)
PyObject *to_python<signed long>::convert(signed long l)
{
  return PyInt_FromLong(l);
}
PYTHONIC_INT_TO_PYTHON(unsigned long long)
PYTHONIC_INT_TO_PYTHON(signed long long)

#undef PYTHONIC_INT_TO_PYTHON

#define PYTHONIC_INT_FROM_PYTHON(TYPE, NTYPE)                                                      \
  inline bool from_python<TYPE>::is_convertible(PyObject *obj)                                     \
  {                                                                                                \
    return PyInt_CheckExact(obj) || PyObject_TypeCheck(obj, &Py##NTYPE##ArrType_Type);             \
  }                                                                                                \
  inline TYPE from_python<TYPE>::convert(PyObject *obj)                                            \
  {                                                                                                \
    return PyInt_AsLong(obj);                                                                      \
  }

PYTHONIC_INT_FROM_PYTHON(unsigned char, UByte)
PYTHONIC_INT_FROM_PYTHON(signed char, Byte)
PYTHONIC_INT_FROM_PYTHON(unsigned short, UShort)
PYTHONIC_INT_FROM_PYTHON(signed short, Short)
PYTHONIC_INT_FROM_PYTHON(unsigned int, UInt)
PYTHONIC_INT_FROM_PYTHON(signed int, Int)
PYTHONIC_INT_FROM_PYTHON(unsigned long, ULong)
PYTHONIC_INT_FROM_PYTHON(signed long, Long)
PYTHONIC_INT_FROM_PYTHON(unsigned long long, ULongLong)
PYTHONIC_INT_FROM_PYTHON(signed long long, LongLong)

#undef PYTHONIC_INT_FROM_PYTHON

PYTHONIC_NS_END
#endif

#endif
