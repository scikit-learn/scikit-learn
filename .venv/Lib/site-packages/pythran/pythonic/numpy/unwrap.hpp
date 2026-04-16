#ifndef PYTHONIC_NUMPY_UNWRAP_HPP
#define PYTHONIC_NUMPY_UNWRAP_HPP

#include "pythonic/include/numpy/unwrap.hpp"

#include "pythonic/numpy/pi.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/int_.hpp"

#include <pythonic/numpy/abs.hpp>
#include <pythonic/numpy/maximum.hpp>
#include <pythonic/numpy/round.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace
  {
    template <class I0, class I1>
    void _unwrap(I0 ibegin, I0 iend, I1 obegin, double discont, utils::int_<1>)
    {
      *obegin = *ibegin;
      ++ibegin;
      for (; ibegin != iend; ++ibegin, ++obegin) {
        if (functor::abs{}(*obegin - *ibegin) > discont)
          *(obegin + 1) = *ibegin + 2 * pi * functor::round{}((*obegin - *ibegin) / (2 * pi));
        else
          *(obegin + 1) = *ibegin;
      }
    }

    template <class I0, class I1, size_t N>
    void _unwrap(I0 ibegin, I0 iend, I1 obegin, double discont, utils::int_<N>)
    {
      for (; ibegin != iend; ++ibegin, ++obegin)
        _unwrap((*ibegin).begin(), (*ibegin).end(), (*obegin).begin(), discont,
                utils::int_<N - 1>());
    }
  } // namespace

  template <class E>
  types::ndarray<double, typename E::shape_t> unwrap(E const &expr, double discont)
  {
    discont = functor::maximum{}(discont, pi);
    types::ndarray<double, typename E::shape_t> out(sutils::getshape(expr), builtins::None);
    _unwrap(expr.begin(), expr.end(), out.begin(), discont, utils::int_<E::value>());
    return out;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
