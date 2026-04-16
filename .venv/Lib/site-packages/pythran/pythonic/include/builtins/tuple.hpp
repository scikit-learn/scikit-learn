#ifndef PYTHONIC_INCLUDE_BUILTIN_TUPLE_HPP
#define PYTHONIC_INCLUDE_BUILTIN_TUPLE_HPP

#include "pythonic/include/types/dynamic_tuple.hpp"
#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  template <class... Types>
  std::tuple<Types...> tuple(std::tuple<Types...> const &t);

  inline std::tuple<> tuple()
  {
    return {};
  }

  template <class Iterable> /* this is far from perfect, but how to cope with
                               the
                               difference between python tuples && c++ ones ? */
      std::enable_if_t <
      types::len_of<std::remove_cv_t<std::remove_reference_t<Iterable>>>::value<
          0, types::dynamic_tuple<typename std::iterator_traits<typename std::remove_cv_t<
                 std::remove_reference_t<Iterable>>::iterator>::value_type>>
      tuple(Iterable &&i);

  template <class StaticIterable> /* specialization if we are capable to statically
                                     compute the size of the input */
  std::enable_if_t<
      types::len_of<std::remove_cv_t<std::remove_reference_t<StaticIterable>>>::value >= 0,
      types::array_tuple<
          typename std::iterator_traits<typename std::remove_cv_t<
              std::remove_reference_t<StaticIterable>>::iterator>::value_type,
          types::len_of<std::remove_cv_t<std::remove_reference_t<StaticIterable>>>::value>>
  tuple(StaticIterable &&i);

  DEFINE_FUNCTOR(pythonic::builtins, tuple);
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<builtins::functor::tuple> {
  static PyObject *convert(builtins::functor::tuple const &c);
};

template <>
struct from_python<builtins::functor::tuple> {
  static bool is_convertible(PyObject *obj);
  static builtins::functor::tuple convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif
#endif
