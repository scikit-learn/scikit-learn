#ifndef PYTHONIC_BUILTIN_BOOL_HPP
#define PYTHONIC_BUILTIN_BOOL_HPP

#include "pythonic/include/builtins/bool_.hpp"

#include "pythonic/types/tuple.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{
  namespace functor
  {

    template <class T>
    bool bool_::operator()(T const &val) const
    {
      return static_cast<bool>(val);
    }

    template <class... Ts>
    bool bool_::operator()(std::tuple<Ts...> const &val) const
    {
      return sizeof...(Ts);
    }

    template <class T, size_t N>
    bool bool_::operator()(types::array_tuple<T, N> const &val) const
    {
      return N;
    }

    inline bool bool_::operator()() const
    {
      return false;
    }
  } // namespace functor
} // namespace builtins
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

inline PyObject *to_python<builtins::functor::bool_>::convert(builtins::functor::bool_ const &c)
{
  return (PyObject *)&PyBool_Type;
}

inline bool from_python<builtins::functor::bool_>::is_convertible(PyObject *obj)
{
  return obj == (PyObject *)&PyBool_Type;
}

inline builtins::functor::bool_ from_python<builtins::functor::bool_>::convert(PyObject *obj)
{
  return {};
}

PYTHONIC_NS_END
#endif

#endif
