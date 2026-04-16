#ifndef PYTHONIC_NUMPY_INTERSECT1D_HPP
#define PYTHONIC_NUMPY_INTERSECT1D_HPP

#include "pythonic/include/numpy/intersect1d.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/combined.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/pdqsort.hpp"

#include <algorithm>
#include <set>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class F>
  types::ndarray<typename __combined<typename E::dtype, typename F::dtype>::type,
                 types::pshape<long>>
  intersect1d(E const &e, F const &f)
  {
    using T = typename __combined<typename E::dtype, typename F::dtype>::type;
    auto ae = asarray(e);
    auto af = asarray(f);
    std::set<T, std::less<T>, utils::allocator<T>> sae(ae.fbegin(), ae.fend());
    std::set<T, std::less<T>, utils::allocator<T>> found;
    types::list<T> lout(0);
    lout.reserve(sae.size());
    for (auto iter = af.fbegin(), end = af.fend(); iter != end; ++iter) {
      auto curr = *iter;
      if (sae.find(curr) != sae.end() && found.find(curr) == found.end()) {
        found.insert(curr);
        lout.push_back(curr);
      }
    }
    pdqsort(lout.begin(), lout.end());
    return {lout};
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
