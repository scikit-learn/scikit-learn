#ifndef PYTHONIC_NUMPY_NANMIN_HPP
#define PYTHONIC_NUMPY_NANMIN_HPP

#include "pythonic/include/numpy/nanmin.hpp"

#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/numpy/isnan.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace
  {
    template <class E, class F>
    bool _nanmin(E begin, E end, F &min, utils::int_<1>)
    {
      bool found = false;
      for (; begin != end; ++begin) {
        auto curr = *begin;
        if (!functor::isnan()(curr) && curr <= min) {
          min = curr;
          found = true;
        }
      }
      return found;
    }

    template <class E, class F, size_t N>
    bool _nanmin(E begin, E end, F &min, utils::int_<N>)
    {
      bool found = false;
      for (; begin != end; ++begin)
        found |= _nanmin((*begin).begin(), (*begin).end(), min, utils::int_<N - 1>());
      return found;
    }
  } // namespace

  template <class E>
  typename E::dtype nanmin(E const &expr)
  {
    bool found = false;
    typename E::dtype min = std::numeric_limits<typename E::dtype>::max();
    found = _nanmin(expr.begin(), expr.end(), min, utils::int_<E::value>());
    if (!found)
      min = std::numeric_limits<typename E::dtype>::quiet_NaN();
    return min;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
