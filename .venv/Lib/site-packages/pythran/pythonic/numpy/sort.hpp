#ifndef PYTHONIC_NUMPY_SORT_HPP
#define PYTHONIC_NUMPY_SORT_HPP

#include "pythonic/include/numpy/sort.hpp"
#include "pythonic/numpy/ndarray/sort.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{

  template <class E>
  types::ndarray<typename E::dtype, types::array_tuple<long, E::value>> sort(E const &expr,
                                                                             long axis)
  {
    auto out = functor::array{}(expr);
    ndarray::sort(out, axis);
    return out;
  }

  template <class E>
  types::ndarray<typename E::dtype, types::array_tuple<long, 1>> sort(E const &expr,
                                                                      types::none_type)
  {
    auto out = functor::array{}(expr).flat();
    ndarray::sort(out, types::none_type{});
    return out;
  }

  template <class E>
  types::ndarray<typename E::dtype, types::array_tuple<long, E::value>>
  sort(E const &expr, long axis, types::str const &kind)
  {
    auto out = functor::array{}(expr);
    ndarray::sort(out, axis, kind);
    return out;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
