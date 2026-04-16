#ifndef PYTHONIC_INCLUDE_BUILTIN_BOOL_HPP
#define PYTHONIC_INCLUDE_BUILTIN_BOOL_HPP

#include "pythonic/include/types/tuple.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace functor
  {

    struct bool_ {
      using callable = void;
      using type = bool;

      bool operator()() const;

      template <class T>
      bool operator()(T const &val) const;

      template <class... Ts>
      bool operator()(std::tuple<Ts...> const &val) const;

      template <class T, size_t N>
      bool operator()(types::array_tuple<T, N> const &val) const;

      friend std::ostream &operator<<(std::ostream &os, bool_)
      {
        return os << "bool";
      }
    };
  } // namespace functor
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <>
struct to_python<builtins::functor::bool_> {
  static PyObject *convert(builtins::functor::bool_ const &c);
};

template <>
struct from_python<builtins::functor::bool_> {
  static bool is_convertible(PyObject *obj);
  static builtins::functor::bool_ convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif
#endif
