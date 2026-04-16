#ifndef PYTHONIC_BUILTIN_SLICE_HPP
#define PYTHONIC_BUILTIN_SLICE_HPP

#include "pythonic/include/builtins/slice.hpp"
#include "pythonic/types/slice.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {
    inline types::cstride_slice<1> slice(types::none<long> stop)
    {
      return {types::none<long>(), stop};
    }
    inline types::cstride_slice<1> slice(types::none<long> start, types::none<long> stop)
    {
      return {start, stop};
    }
    inline types::slice slice(types::none<long> start, types::none<long> stop,
                              types::none<long> step)
    {
      return {start, stop, step};
    }
  } // namespace anonymous
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<builtins::functor::slice>::convert(builtins::functor::slice const &c)
{
  return (PyObject *)&PySlice_Type;
}

inline bool from_python<builtins::functor::slice>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PySlice_Type;
}

inline builtins::functor::slice from_python<builtins::functor::slice>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
