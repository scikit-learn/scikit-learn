#ifndef PYTHONIC_INCLUDE_BUILTIN_SLICE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SLICE_HPP

#include "pythonic/include/types/slice.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {
    types::cstride_slice<1> slice(types::none<long> stop);
    types::cstride_slice<1> slice(types::none<long> start, types::none<long> stop);
    types::slice slice(types::none<long> start, types::none<long> stop, types::none<long> step);
  } // namespace anonymous

  DEFINE_FUNCTOR(pythonic::builtins::anonymous, slice);
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<builtins::functor::slice> {
  static PyObject *convert(builtins::functor::slice const &c);
};

template <>
struct from_python<builtins::functor::slice> {
  static bool is_convertible(PyObject *obj);
  static builtins::functor::slice convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
