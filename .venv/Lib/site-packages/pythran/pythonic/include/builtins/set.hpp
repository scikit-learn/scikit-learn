#ifndef PYTHONIC_INCLUDE_BUILTIN_SET_HPP
#define PYTHONIC_INCLUDE_BUILTIN_SET_HPP

#include "pythonic/include/types/set.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {
    inline types::empty_set set();

    template <class Iterable>
    inline types::set<typename std::iterator_traits<
        typename std::remove_reference_t<Iterable>::iterator>::value_type>
    set(Iterable &&t);
  } // namespace anonymous

  DEFINE_FUNCTOR(pythonic::builtins::anonymous, set);
} // namespace builtins
PYTHONIC_NS_END
#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<builtins::functor::set> {
  static PyObject *convert(builtins::functor::set const &c);
};

template <>
struct from_python<builtins::functor::set> {
  static bool is_convertible(PyObject *obj);
  static builtins::functor::set convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif
#endif
