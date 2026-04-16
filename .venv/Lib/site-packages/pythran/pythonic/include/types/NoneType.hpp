#ifndef PYTHONIC_INCLUDE_TYPES_NONE_HPP
#define PYTHONIC_INCLUDE_TYPES_NONE_HPP

#include "pythonic/include/operator_/mod.hpp"
#include "pythonic/include/types/assignable.hpp"
#include <ostream>

PYTHONIC_NS_BEGIN

namespace types
{

  static const intptr_t NONE_ID = 0x1331;

  struct none_type {
    none_type();
    intptr_t id() const;
    bool operator==(none_type) const
    {
      return true;
    }
    explicit operator bool() const
    {
      return false;
    }
  };

  inline std::ostream &operator<<(std::ostream &os, none_type const &)
  {
    return os << "None";
  }

  template <class T, bool is_fundamental = std::is_fundamental<T>::value>
  struct none;

  /* Type adaptor to simulate an option type
   *
   * see http://en.wikipedia.org/wiki/Option_type
   */
  template <class T>
  struct none<T, false> : T {

    bool is_none; // set to true if the type is none

    none(none_type const &);

    none() : T(), is_none{true}
    {
    }
    none(none const &other) = default;
    none(T const &arg) : T(arg), is_none(false)
    {
    }
    template <class OT>
    none(OT const &arg) : none(T(arg))
    {
    }

    bool operator==(none_type const &) const;

    template <class O>
    bool operator==(O const &t) const;

    bool operator!=(none_type const &) const;

    template <class O>
    bool operator!=(O const &t) const;

    explicit operator bool() const;

    explicit operator T const &() const
    {
      assert(!is_none);
      return *static_cast<T const *>(this);
    }

    intptr_t id() const;
  };

  /* specialization of none for integral types we cannot derive from
   */
  template <class P, class T>
  struct none_data {
    explicit operator bool() const
    {
      return !static_cast<P const *>(this)->is_none && static_cast<P const *>(this)->data;
    }
    operator T const &() const
    {
      return static_cast<P const *>(this)->data;
    }
  };
  template <class P>
  struct none_data<P, bool> {
    operator bool() const
    {
      return !static_cast<P const *>(this)->is_none && static_cast<P const *>(this)->data;
    }
  };
  template <class T>
  struct none<T, true> : none_data<none<T, true>, T> {
    T data;

    template <class T1>
    none &operator+=(T1 other);
    template <class T1>
    none &operator-=(T1 other);
    template <class T1>
    none &operator*=(T1 other);
    template <class T1>
    none &operator/=(T1 other);

    bool is_none;
    none();
    none(none_type const &);
    none(T const &data);
    bool operator==(none_type const &) const;
    template <class O>
    bool operator==(O const &t) const;
    bool operator!=(none_type const &) const;
    template <class O>
    bool operator!=(O const &t) const;
    T &operator=(T const &t);
    intptr_t id() const;
    template <class T1>
    operator none<T1, true>()
    {
      if (is_none)
        return {none_type{}};
      else
        return {static_cast<T1>(data)};
    }
  };

#define NONE_OPERATOR_OVERLOAD(op)                                                                 \
  template <class T>                                                                               \
  auto operator op(none<T> const &t0, T const &t1)->decltype(static_cast<T const &>(t0) op t1)     \
  {                                                                                                \
    return static_cast<T const &>(t0) op t1;                                                       \
  }                                                                                                \
                                                                                                   \
  template <class T>                                                                               \
  auto operator op(T const &t0, none<T> const &t1)->decltype(t0 op static_cast<T const &>(t1))     \
  {                                                                                                \
    return t0 op static_cast<T const &>(t1);                                                       \
  }                                                                                                \
                                                                                                   \
  template <class T>                                                                               \
  auto operator op(none<T> const &t0, none<T> const &t1)                                           \
      ->none<decltype(static_cast<T const &>(t0) op static_cast<T const &>(t1))>                   \
  {                                                                                                \
    if (t0.is_none && t1.is_none)                                                                  \
      return none_type{};                                                                          \
    else {                                                                                         \
      return {static_cast<T const &>(t0) op static_cast<T const &>(t1)};                           \
    }                                                                                              \
  }

  NONE_OPERATOR_OVERLOAD(+)
  NONE_OPERATOR_OVERLOAD(-)
  NONE_OPERATOR_OVERLOAD(*)
  NONE_OPERATOR_OVERLOAD(/)

  NONE_OPERATOR_OVERLOAD(>)
  NONE_OPERATOR_OVERLOAD(>=)
  NONE_OPERATOR_OVERLOAD(<)
  NONE_OPERATOR_OVERLOAD(<=)

  template <class T0, class T1>
  decltype(operator_::mod(std::declval<T0>(), std::declval<T1>())) operator%(none<T0> const &t0,
                                                                             T1 const &t1);

  template <class T0, class T1>
  decltype(operator_::mod(std::declval<T0>(), std::declval<T1>())) operator%(T0 const &t0,
                                                                             none<T1> const &t1);

  template <class T0, class T1>
  none<decltype(operator_::mod(std::declval<T0>(), std::declval<T1>()))>
  operator%(none<T0> const &t0, none<T1> const &t1);

  template <class T, bool F>
  std::ostream &operator<<(std::ostream &os, none<T, F> const &v)
  {
    if (v.is_none)
      return os << none_type();
    else
      return os << v.data;
  }

  template <class T>
  struct is_none {
    static const bool value = false;
  };

  template <class T>
  struct is_none<none<T>> {
    static const bool value = true;
  };
} // namespace types

template <class T>
struct assignable<types::none<T>> {
  using type = types::none<typename assignable<T>::type>;
};
PYTHONIC_NS_END

namespace std
{
  /* std::get overload */
  template <size_t I, class T0>
  auto get(pythonic::types::none<T0> const &t) -> decltype(std::get<I>((T0 const &)t));

  template <size_t I, class T0>
  struct tuple_element<I, pythonic::types::none<T0>> {
    using type = std::tuple_element_t<I, T0>;
  };

  template <>
  struct hash<pythonic::types::none_type> {
    size_t operator()(const pythonic::types::none_type &x) const
    {
      return 0;
    }
  };
} // namespace std

/* type inference stuff { */
#include "pythonic/include/types/combined.hpp"

template <class T0, class T1>
struct __combined<pythonic::types::none<T0>, T1> {
  static_assert(!pythonic::types::is_none<T1>::value, "none of none should'nt exist");
  using type = pythonic::types::none<typename __combined<T0, T1>::type>;
};

template <class T0, class T1>
struct __combined<T1, pythonic::types::none<T0>> {
  static_assert(!pythonic::types::is_none<T0>::value, "none of none should'nt exist");
  using type = pythonic::types::none<typename __combined<T0, T1>::type>;
};

template <class T0, class T1>
struct __combined<pythonic::types::none<T1>, pythonic::types::none<T0>> {
  static_assert(!pythonic::types::is_none<T0>::value, "none of none shouldn't exist");
  static_assert(!pythonic::types::is_none<T1>::value, "none of none shouldn't exist");
  using type = pythonic::types::none<typename __combined<T0, T1>::type>;
};

template <class T>
struct __combined<pythonic::types::none_type, T> {
  static_assert(!pythonic::types::is_none<T>::value, "none of none shouldn't exist");
  using type = pythonic::types::none<T>;
};

template <class T>
struct __combined<pythonic::types::none_type, pythonic::types::none<T>> {
  static_assert(!pythonic::types::is_none<T>::value, "none of none shouldn't exist");
  using type = pythonic::types::none<T>;
};

template <class T>
struct __combined<T, pythonic::types::none_type> {
  static_assert(!pythonic::types::is_none<T>::value, "none of none shouldn't exist");
  using type = pythonic::types::none<T>;
};
template <class T>
struct __combined<pythonic::types::none<T>, pythonic::types::none_type> {
  static_assert(!pythonic::types::is_none<T>::value, "none of none shouldn't exist");
  using type = pythonic::types::none<T>;
};

template <>
struct __combined<pythonic::types::none_type, pythonic::types::none_type> {
  using type = pythonic::types::none_type;
};

/* } */

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN
template <>
struct to_python<types::none_type> {
  static PyObject *convert(types::none_type);
};

template <class T>
struct to_python<types::none<T>> {
  static PyObject *convert(types::none<T> const &n);
};

template <>
struct from_python<types::none_type> {

  static bool is_convertible(PyObject *obj);

  static types::none_type convert(PyObject *obj);
};
PYTHONIC_NS_END

#endif
#endif
