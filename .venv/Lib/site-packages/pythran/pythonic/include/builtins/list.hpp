#ifndef PYTHONIC_INCLUDE_BUILTIN_LIST_HPP
#define PYTHONIC_INCLUDE_BUILTIN_LIST_HPP

#include "pythonic/include/types/list.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <iterator>
#include <type_traits>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {

    inline types::empty_list list();
    inline types::empty_list list(types::empty_list);

    template <class Iterable>
    types::list<std::decay_t<typename std::iterator_traits<
        typename std::remove_reference_t<Iterable>::iterator>::value_type>>
    list(Iterable &&t);
  } // namespace anonymous

  DEFINE_FUNCTOR(pythonic::builtins::anonymous, list);
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<builtins::functor::list> {
  static PyObject *convert(builtins::functor::list const &c);
};

template <>
struct from_python<builtins::functor::list> {
  static bool is_convertible(PyObject *obj);
  static builtins::functor::list convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
