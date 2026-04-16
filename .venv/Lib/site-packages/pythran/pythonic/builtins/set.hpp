#ifndef PYTHONIC_BUILTIN_SET_HPP
#define PYTHONIC_BUILTIN_SET_HPP

#include "pythonic/include/builtins/set.hpp"

#include "pythonic/types/set.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {

    inline types::empty_set set()
    {
      return types::empty_set();
    }

    template <class Iterable>
    inline types::set<typename std::iterator_traits<
        typename std::remove_reference_t<Iterable>::iterator>::value_type>
    set(Iterable &&t)
    {
      return {t.begin(), t.end()};
    }
  } // namespace anonymous
} // namespace builtins
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<builtins::functor::set>::convert(builtins::functor::set const &c)
{
  return (PyObject *)&PySet_Type;
}

inline bool from_python<builtins::functor::set>::is_convertible(PyObject *obj)
{
  if (obj == (PyObject *)&PySet_Type)
    return true;
  PyObject *Origin = PyObject_GetAttrString(obj, "__origin__");
  if (!Origin)
    return false;
  bool res = (Origin == (PyObject *)&PySet_Type);
  Py_DECREF(Origin);
  return res;
}

inline builtins::functor::set from_python<builtins::functor::set>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif
#endif
