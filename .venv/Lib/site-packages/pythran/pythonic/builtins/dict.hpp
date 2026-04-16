#ifndef PYTHONIC_BUILTIN_DICT_HPP
#define PYTHONIC_BUILTIN_DICT_HPP

#include "pythonic/include/builtins/dict.hpp"

#include "pythonic/types/dict.hpp"
#include "pythonic/utils/functor.hpp"

#include <tuple>
#include <utility>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {
    inline types::empty_dict dict()
    {
      return types::empty_dict();
    }

    template <class K, class V>
    types::dict<K, V> dict(types::dict<K, V> const &other)
    {
      return other.copy();
    }

    template <class Iterable>
    auto dict(Iterable &&iterable)
        -> types::dict<std::decay_t<decltype(std::get<0>(*iterable.begin()))>,
                       std::decay_t<decltype(std::get<1>(*iterable.begin()))>>
    {
      types::dict<std::decay_t<decltype(std::get<0>(*iterable.begin()))>,
                  std::decay_t<decltype(std::get<1>(*iterable.begin()))>>
          out = types::empty_dict();
      for (auto const &i : iterable)
        out[std::get<0>(i)] = std::get<1>(i);
      return out;
    }
  } // namespace anonymous
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<builtins::functor::dict>::convert(builtins::functor::dict const &c)
{
  return (PyObject *)&PyDict_Type;
}

inline bool from_python<builtins::functor::dict>::is_convertible(PyObject *obj)
{
  if (obj == (PyObject *)&PyDict_Type)
    return true;
  PyObject *Origin = PyObject_GetAttrString(obj, "__origin__");
  if (!Origin)
    return false;
  bool res = (Origin == (PyObject *)&PyDict_Type);
  Py_DECREF(Origin);
  return res;
}

inline builtins::functor::dict from_python<builtins::functor::dict>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
