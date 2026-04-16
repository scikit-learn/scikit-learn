#ifndef PYTHONIC_INCLUDE_UTILS_ARRAY_HELPER_HPP
#define PYTHONIC_INCLUDE_UTILS_ARRAY_HELPER_HPP

#include "pythonic/include/types/tuple.hpp"

PYTHONIC_NS_BEGIN

/* recursively return the value at the position given by `indices' in
 * the `self' "array like". It may be a sub array instead of real value.
 * indices[0] is the coordinate for the first dimension && indices[M-1]
 * is for the last one.
 */
template <size_t L>
struct nget {
  template <class A, size_t M>
  auto operator()(A &&self, types::array_tuple<long, M> const &indices)
      -> decltype(nget<L - 1>()(std::forward<A>(self)[0], indices));
  template <class A, size_t M>
  auto fast(A &&self, types::array_tuple<long, M> const &indices)
      -> decltype(nget<L - 1>().fast(std::forward<A>(self).fast(0), indices));
};

template <>
struct nget<0> {
  template <class A, size_t M>
  auto operator()(A &&self, types::array_tuple<long, M> const &indices)
      -> decltype(std::forward<A>(self)[indices[M - 1]]);
  template <class A, size_t M>
  auto fast(A &&self, types::array_tuple<long, M> const &indices)
      -> decltype(std::forward<A>(self).fast(indices[M - 1]));
};
PYTHONIC_NS_END

#endif
