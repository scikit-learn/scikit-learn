#ifndef PYTHONIC_BUILTIN_LIST_HPP
#define PYTHONIC_BUILTIN_LIST_HPP

#include "pythonic/include/builtins/list.hpp"

#include "pythonic/types/list.hpp"
#include "pythonic/utils/functor.hpp"

#include <iterator>
#include <type_traits>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {

    inline types::empty_list list()
    {
      return types::empty_list();
    }

    inline types::empty_list list(types::empty_list)
    {
      return types::empty_list();
    }

    template <class Iterable>
    types::list<std::decay_t<typename std::iterator_traits<
        typename std::remove_reference_t<Iterable>::iterator>::value_type>>
    list(Iterable &&t)
    {
      return types::list<std::decay_t<typename std::iterator_traits<
          typename std::remove_reference_t<Iterable>::iterator>::value_type>>(t.begin(), t.end());
    }
  } // namespace anonymous
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<builtins::functor::list>::convert(builtins::functor::list const &c)
{
  return (PyObject *)&PyList_Type;
}

inline bool from_python<builtins::functor::list>::is_convertible(PyObject *obj)
{
  if (obj == (PyObject *)&PyList_Type)
    return true;
  PyObject *Origin = PyObject_GetAttrString(obj, "__origin__");
  if (!Origin)
    return false;
  bool res = (Origin == (PyObject *)&PyList_Type);
  Py_DECREF(Origin);
  return res;
}

inline builtins::functor::list from_python<builtins::functor::list>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
