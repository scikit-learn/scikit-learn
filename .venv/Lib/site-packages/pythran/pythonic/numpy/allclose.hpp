#ifndef PYTHONIC_NUMPY_ALLCLOSE_HPP
#define PYTHONIC_NUMPY_ALLCLOSE_HPP

#include "pythonic/include/numpy/allclose.hpp"

#include "pythonic/numpy/abs.hpp"
#include "pythonic/numpy/isfinite.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace
  {
    template <class I0, class I1>
    bool _allclose(I0 begin, I0 end, I1 ibegin, double rtol, double atol, utils::int_<1>)
    {
      for (; begin != end; ++begin, ++ibegin) {
        auto u = *begin;
        auto v = *ibegin;
        if (((!functor::isfinite()(u) || !functor::isfinite()(v)) &&
             u != v) || // Infinite && NaN cases
            functor::abs()(u - v) > (atol + rtol * functor::abs()(v))) {
          return false;
        }
      }
      return true;
    }

    template <class I0, class I1, size_t N>
    bool _allclose(I0 begin, I0 end, I1 ibegin, double rtol, double atol, utils::int_<N>)
    {
      for (; begin != end; ++begin, ++ibegin)
        if (!_allclose((*begin).begin(), (*begin).end(), (*ibegin).begin(), rtol, atol,
                       utils::int_<N - 1>()))
          return false;
      return true;
    }
  } // namespace

  template <class U, class V>
  bool allclose(U const &u, V const &v, double rtol, double atol)
  {
    return _allclose(u.begin(), u.end(), v.begin(), rtol, atol, utils::int_<U::value>());
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
