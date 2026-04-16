#ifndef PYTHONIC_UTILS_ARRAY_HELPER_HPP
#define PYTHONIC_UTILS_ARRAY_HELPER_HPP

#include "pythonic/include/utils/array_helper.hpp"
#include "pythonic/types/tuple.hpp"

PYTHONIC_NS_BEGIN

/* recursively return the value at the position given by `indices' in the
 * `self' "array like". It may be a sub array instead of real value.
 * indices[0] is the coordinate for the first dimension && indices[M-1] is
 * for the last one.
 */
template <size_t L>
template <class A, size_t M>
auto nget<L>::operator()(A &&self, types::array_tuple<long, M> const &indices)
    -> decltype(nget<L - 1>()(std::forward<A>(self)[0], indices))
{
  return nget<L - 1>()(std::forward<A>(self)[indices[M - L - 1]], indices);
}

template <size_t L>
template <class A, size_t M>
auto nget<L>::fast(A &&self, types::array_tuple<long, M> const &indices)
    -> decltype(nget<L - 1>().fast(std::forward<A>(self).fast(0), indices))
{
  return nget<L - 1>().fast(std::forward<A>(self).fast(indices[M - L - 1]), indices);
}

template <class A, size_t M>
auto nget<0>::operator()(A &&self, types::array_tuple<long, M> const &indices)
    -> decltype(std::forward<A>(self)[indices[M - 1]])
{
  return std::forward<A>(self)[indices[M - 1]];
}

template <class A, size_t M>
auto nget<0>::fast(A &&self, types::array_tuple<long, M> const &indices)
    -> decltype(std::forward<A>(self).fast(indices[M - 1]))
{
  return std::forward<A>(self).fast(indices[M - 1]);
}
PYTHONIC_NS_END
#endif
