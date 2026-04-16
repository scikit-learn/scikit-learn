#ifndef PYTHONIC_BUILTIN_TUPLE_HPP
#define PYTHONIC_BUILTIN_TUPLE_HPP

#include "pythonic/include/builtins/tuple.hpp"

#include "pythonic/types/dynamic_tuple.hpp"
#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class... Types>
  std::tuple<Types...> tuple(std::tuple<Types...> const &t)
  {
    return t;
  }

  template <class Iterable>
      /* this is far from perfect, but how to cope with the
         difference between python tuples and c++ ones ? */
      std::enable_if_t <
      types::len_of<std::remove_cv_t<std::remove_reference_t<Iterable>>>::value<
          0, types::dynamic_tuple<typename std::iterator_traits<typename std::remove_cv_t<
                 std::remove_reference_t<Iterable>>::iterator>::value_type>>
      tuple(Iterable &&i)
  {
    return {i.begin(), i.end()};
  }

  template <class StaticIterable>
  /* specialization if we are capable to statically compute the size of the
     input */
  std::enable_if_t<
      types::len_of<std::remove_cv_t<std::remove_reference_t<StaticIterable>>>::value >= 0,
      types::array_tuple<
          typename std::iterator_traits<typename std::remove_cv_t<
              std::remove_reference_t<StaticIterable>>::iterator>::value_type,
          types::len_of<std::remove_cv_t<std::remove_reference_t<StaticIterable>>>::value>>
  tuple(StaticIterable &&i)
  {
    types::array_tuple<
        typename std::iterator_traits<typename std::remove_cv_t<
            std::remove_reference_t<StaticIterable>>::iterator>::value_type,
        types::len_of<std::remove_cv_t<std::remove_reference_t<StaticIterable>>>::value>
        res;
    std::copy(i.begin(), i.end(), res.begin());
    return res;
  }
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<builtins::functor::tuple>::convert(builtins::functor::tuple const &c)
{
  return (PyObject *)&PyTuple_Type;
}

inline bool from_python<builtins::functor::tuple>::is_convertible(PyObject *obj)
{
  if (obj == (PyObject *)&PyTuple_Type)
    return true;
  PyObject *Origin = PyObject_GetAttrString(obj, "__origin__");
  if (!Origin)
    return false;
  bool res = (Origin == (PyObject *)&PyTuple_Type);
  Py_DECREF(Origin);
  return res;
}

inline builtins::functor::tuple from_python<builtins::functor::tuple>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
