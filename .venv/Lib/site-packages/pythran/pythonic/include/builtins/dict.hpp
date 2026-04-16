#ifndef PYTHONIC_INCLUDE_BUILTIN_DICT_HPP
#define PYTHONIC_INCLUDE_BUILTIN_DICT_HPP

#include "pythonic/include/types/dict.hpp"
#include "pythonic/include/utils/functor.hpp"

#include <tuple>
#include <utility>

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {
    types::empty_dict dict();

    template <class K, class V>
    types::dict<K, V> dict(types::dict<K, V> const &);

    template <class Iterable>
    auto dict(Iterable &&iterable)
        -> types::dict<std::decay_t<decltype(std::get<0>(*iterable.begin()))>,
                       std::decay_t<decltype(std::get<1>(*iterable.begin()))>>;
  } // namespace anonymous

  DEFINE_FUNCTOR(pythonic::builtins::anonymous, dict);
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<builtins::functor::dict> {
  static PyObject *convert(builtins::functor::dict const &c);
};

template <>
struct from_python<builtins::functor::dict> {
  static bool is_convertible(PyObject *obj);
  static builtins::functor::dict convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif
#endif
