#ifndef PYTHONIC_INCLUDE_BUILTIN_FLOAT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_FLOAT_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace functor
  {
    struct float_ {
      using callable = void;
      using type = double;

      template <class T>
      type operator()(T &&t) const;

      type operator()() const;

      friend std::ostream &operator<<(std::ostream &os, float_)
      {
        return os << "float";
      }
    };
  } // namespace functor
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<builtins::functor::float_> {
  static PyObject *convert(builtins::functor::float_ const &c);
};

template <>
struct from_python<builtins::functor::float_> {
  static bool is_convertible(PyObject *obj);
  static builtins::functor::float_ convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
