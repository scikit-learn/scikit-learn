#ifndef PYTHONIC_TYPES_TYPE_HPP
#define PYTHONIC_TYPES_TYPE_HPP

#include "pythonic/include/types/type.hpp"

#include "pythonic/builtins/bool_.hpp"
#include "pythonic/builtins/complex.hpp"
#include "pythonic/builtins/dict.hpp"
#include "pythonic/builtins/float_.hpp"
#include "pythonic/builtins/int_.hpp"
#include "pythonic/builtins/list.hpp"
#include "pythonic/builtins/set.hpp"
#include "pythonic/builtins/slice.hpp"
#include "pythonic/builtins/str.hpp"
#include "pythonic/builtins/tuple.hpp"
#include "pythonic/numpy/array.hpp"
#include "pythonic/numpy/byte.hpp"
#include "pythonic/numpy/complex256.hpp"
#include "pythonic/numpy/complex64.hpp"
#include "pythonic/numpy/float128.hpp"
#include "pythonic/numpy/float32.hpp"
#include "pythonic/numpy/int_.hpp"
#include "pythonic/numpy/intc.hpp"
#include "pythonic/numpy/longlong.hpp"
#include "pythonic/numpy/short_.hpp"
#include "pythonic/numpy/ubyte.hpp"
#include "pythonic/numpy/uint.hpp"
#include "pythonic/numpy/uintc.hpp"
#include "pythonic/numpy/ulonglong.hpp"
#include "pythonic/numpy/ushort.hpp"

PYTHONIC_NS_BEGIN

namespace types
{
  template <>
  struct type_functor<bool> {
    using type = builtins::functor::bool_;
  };
  template <>
  struct type_functor<double> {
    using type = builtins::functor::float_;
  };
  template <>
  struct type_functor<types::str> {
    using type = builtins::functor::str;
  };
  template <>
  struct type_functor<std::complex<float>> {
    using type = numpy::functor::complex64;
  };
  template <>
  struct type_functor<std::complex<double>> {
    using type = builtins::functor::complex;
  };
  template <>
  struct type_functor<std::complex<long double>> {
    using type = numpy::functor::complex256;
  };
  template <>
  struct type_functor<types::empty_set> {
    using type = builtins::functor::set;
  };
  template <class T>
  struct type_functor<types::set<T>> {
    using type = builtins::functor::set;
  };
  template <>
  struct type_functor<types::slice> {
    using type = builtins::functor::slice;
  };
  template <>
  struct type_functor<types::empty_list> {
    using type = builtins::functor::list;
  };
  template <class T>
  struct type_functor<types::list<T>> {
    using type = builtins::functor::list;
  };
  template <class T, size_t N>
  struct type_functor<types::static_list<T, N>> {
    using type = builtins::functor::list;
  };
  template <>
  struct type_functor<types::empty_dict> {
    using type = builtins::functor::dict;
  };
  template <class K, class V>
  struct type_functor<types::dict<K, V>> {
    using type = builtins::functor::dict;
  };
  template <class... Tys>
  struct type_functor<std::tuple<Tys...>> {
    using type = builtins::functor::tuple;
  };
  template <class T, size_t N>
  struct type_functor<types::array_tuple<T, N>> {
    using type = builtins::functor::tuple;
  };
  template <class T, class pS>
  struct type_functor<types::ndarray<T, pS>> {
    using type = numpy::functor::array;
  };
  template <>
  struct type_functor<signed char> {
    using type = numpy::functor::byte;
  };
  template <>
  struct type_functor<unsigned char> {
    using type = numpy::functor::ubyte;
  };
  template <>
  struct type_functor<short> {
    using type = numpy::functor::short_;
  };
  template <>
  struct type_functor<unsigned short> {
    using type = numpy::functor::ushort;
  };
  template <>
  struct type_functor<int> {
    using type = numpy::functor::intc;
  };
  template <>
  struct type_functor<unsigned int> {
    using type = numpy::functor::uintc;
  };
  template <>
  struct type_functor<long> {
    using type = numpy::functor::int_;
  };
  template <>
  struct type_functor<unsigned long> {
    using type = numpy::functor::uint;
  };
  template <>
  struct type_functor<long long> {
    using type = numpy::functor::longlong;
  };
  template <>
  struct type_functor<unsigned long long> {
    using type = numpy::functor::ulonglong;
  };
  template <>
  struct type_functor<float> {
    using type = numpy::functor::float32;
  };
  template <>
  struct type_functor<long double> {
    using type = numpy::functor::float128;
  };

  template <class T>
  typename type_functor<T>::type type(T const &)
  {
    return {};
  }
} // namespace types
PYTHONIC_NS_END

#endif
