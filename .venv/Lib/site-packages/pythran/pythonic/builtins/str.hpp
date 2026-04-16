#ifndef PYTHONIC_BUILTIN_STR_HPP
#define PYTHONIC_BUILTIN_STR_HPP

#include "pythonic/include/builtins/str.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

#include <sstream>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {

    template <class T>
    types::str str(T const &t)
    {
      std::ostringstream oss;
      oss << t;
      return oss.str();
    }

    inline types::str str()
    {
      return "";
    }

    inline types::str str(bool b)
    {
      static char const repr[2][6] = {"False", "True\0"};
      return repr[b];
    }

    inline types::str str(long value)
    {
      /* adapted from http://www.jb.man.ac.uk/~slowe/cpp/itoa.html#performance
       */

      // this buffer is large enough to hold the binary representation, so
      // the decimal representation will be ok
      char buffer[8 * (1 << sizeof(value))];
      char *ptr = buffer, *ptr1 = buffer, tmp_char;
      long tmp_value;

      do {
        tmp_value = value;
        value /= 10;
        *ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmn"
                 "opqrstuvwxyz"[35 + (tmp_value - value * 10)];
      } while (value);

      // Apply negative sign
      if (tmp_value < 0)
        *ptr++ = '-';
      *ptr-- = '\0';
      while (ptr1 < ptr) {
        tmp_char = *ptr;
        *ptr-- = *ptr1;
        *ptr1++ = tmp_char;
      }
      return buffer;
    }

    inline types::str str(double l)
    {
      // when using %g, only 6 significant bits are used, so this should be
      // enough.
      // Use snprintf though
      char buffer[8 * (1 << sizeof(l))];
      snprintf(buffer, sizeof(buffer), "%g", l);
      return buffer;
    }
  } // namespace anonymous
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<builtins::functor::str>::convert(builtins::functor::str const &c)
{
  return (PyObject *)&PyUnicode_Type;
}

inline bool from_python<builtins::functor::str>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyUnicode_Type;
}

inline builtins::functor::str from_python<builtins::functor::str>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
