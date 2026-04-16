#ifndef PYTHONIC_INCLUDE_BUILTIN_STR_HPP
#define PYTHONIC_INCLUDE_BUILTIN_STR_HPP

#include "pythonic/include/types/str.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {
    template <class T>
    types::str str(T const &t);

    inline types::str str();
    inline types::str str(bool b);
    inline types::str str(long value);
    inline types::str str(double l);
  } // namespace anonymous

  DEFINE_FUNCTOR(pythonic::builtins::anonymous, str);
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<builtins::functor::str> {
  static PyObject *convert(builtins::functor::str const &c);
};

template <>
struct from_python<builtins::functor::str> {
  static bool is_convertible(PyObject *obj);
  static builtins::functor::str convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
