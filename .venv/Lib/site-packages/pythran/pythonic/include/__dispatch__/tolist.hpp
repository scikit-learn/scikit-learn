#ifndef PYTHONIC_INCLUDE_DISPATCH_TOLIST_HPP
#define PYTHONIC_INCLUDE_DISPATCH_TOLIST_HPP

#include "pythonic/include/numpy/ndarray/tolist.hpp"
#include "pythonic/include/types/array.hpp"

PYTHONIC_NS_BEGIN

namespace __dispatch__
{
  template <class Any>
  auto tolist(Any &&any) -> decltype(numpy::ndarray::tolist(any))
  {
    return numpy::ndarray::tolist(any);
  }

  template <class T, class S>
  types::list<std::conditional_t<std::is_integral<T>::value, long, double>>
  tolist(types::sliced_array<T, S> &&a)
  {
    return {a.begin(), a.end()};
  }

  template <class T, class S>
  types::list<std::conditional_t<std::is_integral<T>::value, long, double>>
  tolist(types::sliced_array<T, S> const &a)
  {
    return {a.begin(), a.end()};
  }

  template <class T>
  types::list<std::conditional_t<std::is_integral<T>::value, long, double>>
  tolist(types::array<T> &&a)
  {
    return {a.begin(), a.end()};
  }

  template <class T>
  types::list<std::conditional_t<std::is_integral<T>::value, long, double>>
  tolist(types::array<T> const &a)
  {
    return {a.begin(), a.end()};
  }

  DEFINE_FUNCTOR(pythonic::__dispatch__, tolist);
} // namespace __dispatch__
PYTHONIC_NS_END

#endif
