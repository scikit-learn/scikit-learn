#ifndef PYTHONIC_INCLUDE_TYPES_COMPLEX_HPP
#define PYTHONIC_INCLUDE_TYPES_COMPLEX_HPP

#include <complex>

#if defined(_OPENMP)
#pragma omp declare reduction(+ : std::complex<float> : omp_out += omp_in)
#pragma omp declare reduction(* : std::complex<float> : omp_out *= omp_in)
#pragma omp declare reduction(+ : std::complex<double> : omp_out += omp_in)
#pragma omp declare reduction(* : std::complex<double> : omp_out *= omp_in)
#pragma omp declare reduction(+ : std::complex<long double> : omp_out += omp_in)
#pragma omp declare reduction(* : std::complex<long double> : omp_out *= omp_in)
#endif

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace functor
  {
    struct complex64;
    struct complex128;
    struct complex256;
  } // namespace functor
} // namespace numpy

PYTHONIC_NS_END

namespace std
{

  template <class T, class S>
  using complex_broadcast_t =
      std::enable_if_t<std::is_scalar<S>::value && !std::is_same<T, S>::value,
                       std::complex<std::common_type_t<T, S>>>;
  template <class T, class S>
  using complex_bool_t =
      std::enable_if_t<std::is_scalar<S>::value && !std::is_same<T, S>::value, bool>;

  template <class T, class S>
  complex_broadcast_t<T, S> operator+(std::complex<T> self, S other);
  template <class T, class S>
  complex_broadcast_t<T, S> operator+(S self, std::complex<T> other);
  template <class T, class S>
  complex_broadcast_t<T, S> operator-(std::complex<T> self, S other);
  template <class T, class S>
  complex_broadcast_t<T, S> operator-(S self, std::complex<T> other);
  template <class T, class S>
  complex_broadcast_t<T, S> operator*(std::complex<T> self, S other);
  template <class T, class S>
  complex_broadcast_t<T, S> operator*(S self, std::complex<T> other);
  template <class T, class S>
  complex_broadcast_t<T, S> operator/(std::complex<T> self, S other);
  template <class T, class S>
  complex_broadcast_t<T, S> operator/(S self, std::complex<T> other);

  template <class T, class S>
  complex_bool_t<T, S> operator==(std::complex<T> self, S other);
  template <class T, class S>
  complex_bool_t<T, S> operator==(S self, std::complex<T> other);
  template <class T, class S>
  complex_bool_t<T, S> operator!=(std::complex<T> self, S other);
  template <class T, class S>
  complex_bool_t<T, S> operator!=(S self, std::complex<T> other);

  template <class T, class S>
  bool operator<(std::complex<T> self, std::complex<S> other);
  template <class T, class S>
  bool operator<=(std::complex<T> self, std::complex<S> other);
  template <class T, class S>
  bool operator>(std::complex<T> self, std::complex<S> other);
  template <class T, class S>
  bool operator>=(std::complex<T> self, std::complex<S> other);
  template <class T, class S>
  bool operator&&(std::complex<T> self, std::complex<S> other);
  template <class T, class S>
  bool operator||(std::complex<T> self, std::complex<S> other);

  template <class T>
  bool operator!(std::complex<T> self);

  template <class T>
  struct hash<std::complex<T>> {
    size_t operator()(std::complex<T> const &x) const;
  };
} // namespace std

PYTHONIC_NS_BEGIN
namespace builtins
{
  template <class T>
  T getattr(types::attr::REAL, std::complex<T> const &self);
  template <class T>
  T getattr(types::attr::IMAG, std::complex<T> const &self);
  numpy::functor::complex64 getattr(types::attr::DTYPE, std::complex<float> const &self);
  numpy::functor::complex128 getattr(types::attr::DTYPE, std::complex<double> const &self);
  numpy::functor::complex256 getattr(types::attr::DTYPE, std::complex<long double> const &self);
} // namespace builtins
PYTHONIC_NS_END

/* for type inference { */

#include "pythonic/include/types/combined.hpp"
template <class K, class T>
struct __combined<indexable<K>, std::complex<T>> {
  using type = std::complex<T>;
};

template <class K, class T>
struct __combined<std::complex<T>, indexable<K>> {
  using type = std::complex<T>;
};
template <class T0, class T1>
struct __combined<std::complex<T0>, std::complex<T1>> {
  using type = std::complex<typename __combined<T0, T1>::type>;
};

/* } */

#define STD_COMPLEX_IMPLICT_OPERATOR_CAST(op)                                                      \
  template <class T, class U>                                                                      \
  auto operator op(std::complex<T> const &lhs, std::complex<U> const &rhs)                         \
      ->std::complex<std::common_type_t<T, U>>                                                     \
  {                                                                                                \
    using ctype = std::complex<std::common_type_t<T, U>>;                                          \
    return ctype{lhs} op ctype{rhs};                                                               \
  }

STD_COMPLEX_IMPLICT_OPERATOR_CAST(+)
STD_COMPLEX_IMPLICT_OPERATOR_CAST(-)
STD_COMPLEX_IMPLICT_OPERATOR_CAST(*)
STD_COMPLEX_IMPLICT_OPERATOR_CAST(/)
STD_COMPLEX_IMPLICT_OPERATOR_CAST(==)
STD_COMPLEX_IMPLICT_OPERATOR_CAST(!=)

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN

template <class T>
struct to_python<std::complex<T>> {
  static PyObject *convert(std::complex<T> const &c);
};

template <class T>
struct from_python<std::complex<T>> {
  static bool is_convertible(PyObject *obj);
  static std::complex<T> convert(PyObject *obj);
};
PYTHONIC_NS_END
#endif

#endif
