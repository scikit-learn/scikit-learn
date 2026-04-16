#ifndef PYTHONIC_TYPES_NONE_HPP
#define PYTHONIC_TYPES_NONE_HPP

#include "pythonic/include/types/NoneType.hpp"

#include "pythonic/builtins/bool_.hpp"
#include "pythonic/builtins/id.hpp"
#include "pythonic/operator_/mod.hpp"
#include "pythonic/types/assignable.hpp"

PYTHONIC_NS_BEGIN

namespace types
{

  /// none_type implementation
  none_type::none_type()
  {
  }

  inline intptr_t none_type::id() const
  {
    return NONE_ID;
  }

  /// none implementation
  /* Type adapator to simulate an option type
   *
   * see http://en.wikipedia.org/wiki/Option_type
   */
  template <class T>
  none<T, false>::none(none_type const &) : T(), is_none(true)
  {
  }

  template <class T>
  bool none<T, false>::operator==(none_type const &) const
  {
    return is_none;
  }

  template <class T>
  template <class O>
  bool none<T, false>::operator==(O const &t) const
  {
    return !is_none && static_cast<const T &>(*this) == t;
  }

  template <class T>
  bool none<T, false>::operator!=(none_type const &other) const
  {
    return !(*this == other);
  }

  template <class T>
  template <class O>
  bool none<T, false>::operator!=(O const &other) const
  {
    return !(*this == other);
  }

  template <class T>
  none<T, false>::operator bool() const
  {
    return !is_none && builtins::functor::bool_{}(static_cast<const T &>(*this));
  }

  template <class T>
  intptr_t none<T, false>::id() const
  {
    return is_none ? NONE_ID : builtins::id(static_cast<const T &>(*this));
  }
  template <class T>
  std::ostream &operator<<(std::ostream &os, none<T, false> const &v)
  {
    if (v.is_none)
      return os << none_type();
    else
      return os << static_cast<T const &>(v);
  }

  /* specialization of none for integral types we cannot derive from */
  template <class T>
  none<T, true>::none() : data(), is_none(false)
  {
  }

  template <class T>
  none<T, true>::none(none_type const &) : data(), is_none(true)
  {
  }

  template <class T>
  none<T, true>::none(T const &data) : data(data), is_none(false)
  {
  }

  template <class T>
  bool none<T, true>::operator==(none_type const &) const
  {
    return is_none;
  }

  template <class T>
  template <class O>
  bool none<T, true>::operator==(O const &t) const
  {
    return !is_none && data == t;
  }

  template <class T>
  bool none<T, true>::operator!=(none_type const &other) const
  {
    return !(*this == other);
  }

  template <class T>
  template <class O>
  bool none<T, true>::operator!=(O const &other) const
  {
    return !(*this == other);
  }

  template <class T>
  T &none<T, true>::operator=(T const &t)
  {
    is_none = false;
    return data = t;
  }

  template <class T>
  intptr_t none<T, true>::id() const
  {
    return is_none ? NONE_ID : reinterpret_cast<intptr_t>(&data);
  }

  template <class T0, class T1>
  decltype(operator_::mod(std::declval<T0>(), std::declval<T1>())) operator%(none<T0> const &t0,
                                                                             T1 const &t1)
  {
    return operator_::mod(t0.data, t1);
  }

  template <class T0, class T1>
  decltype(operator_::mod(std::declval<T0>(), std::declval<T1>())) operator%(T0 const &t0,
                                                                             none<T1> const &t1)
  {
    return operator_::mod(t0, t1.data);
  }

  template <class T0, class T1>
  none<decltype(operator_::mod(std::declval<T0>(), std::declval<T1>())), true>
  operator%(none<T0> const &t0, none<T1> const &t1)
  {
    if (t0.is_none && t1.is_none)
      return none_type{};
    else
      return {operator_::mod(t0, t1.data)};
  }

  template <class T>
  template <class T1>
  none<T, true> &none<T, true>::operator+=(T1 other)
  {
    if (!is_none)
      data += other;
    return *this;
  }

  template <class T>
  template <class T1>
  none<T, true> &none<T, true>::operator-=(T1 other)
  {
    if (!is_none)
      data -= other;
    return *this;
  }

  template <class T>
  template <class T1>
  none<T, true> &none<T, true>::operator*=(T1 other)
  {
    if (!is_none)
      data *= other;
    return *this;
  }

  template <class T>
  template <class T1>
  none<T, true> &none<T, true>::operator/=(T1 other)
  {
    if (!is_none)
      data /= other;
    return *this;
  }

} // namespace types
PYTHONIC_NS_END

#ifdef ENABLE_PYTHON_MODULE

PYTHONIC_NS_BEGIN

inline bool from_python<types::none_type>::is_convertible(PyObject *obj)
{
  return obj == Py_None;
}

inline types::none_type from_python<types::none_type>::convert(PyObject *obj)
{
  return {};
}

inline PyObject *to_python<types::none_type>::convert(types::none_type)
{
  Py_RETURN_NONE;
}

template <class T>
PyObject *to_python<types::none<T>>::convert(types::none<T> const &n)
{
  if (n.is_none) {
    Py_RETURN_NONE;
  } else {
    return ::to_python(static_cast<T const &>(n));
  }
}
PYTHONIC_NS_END

#endif

#endif
