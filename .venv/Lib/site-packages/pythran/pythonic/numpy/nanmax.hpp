#ifndef PYTHONIC_NUMPY_NANMAX_HPP
#define PYTHONIC_NUMPY_NANMAX_HPP

#include "pythonic/include/numpy/nanmax.hpp"

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
    bool _nanmax(E begin, E end, F &max, utils::int_<1>)
    {
      bool found = false;
      for (; begin != end; ++begin) {
        auto curr = *begin;
        if (!functor::isnan()(curr) && curr >= max) {
          max = curr;
          found = true;
        }
      }
      return found;
    }

    template <class E, class F, size_t N>
    bool _nanmax(E begin, E end, F &max, utils::int_<N>)
    {
      bool found = false;
      for (; begin != end; ++begin)
        found |= _nanmax((*begin).begin(), (*begin).end(), max, utils::int_<N - 1>());
      return found;
    }
  } // namespace

  template <class E>
  typename E::dtype nanmax(E const &expr)
  {
    bool found = false;
    typename E::dtype max = std::numeric_limits<typename E::dtype>::lowest();
    found = _nanmax(expr.begin(), expr.end(), max, utils::int_<E::value>());
    if (!found)
      max = std::numeric_limits<typename E::dtype>::quiet_NaN();
    return max;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
