#ifndef PYTHONIC_NUMPY_SEARCHSORTED_HPP
#define PYTHONIC_NUMPY_SEARCHSORTED_HPP

#include "pythonic/include/numpy/searchsorted.hpp"

#include "pythonic/builtins/None.hpp"
#include "pythonic/builtins/ValueError.hpp"
#include "pythonic/numpy/asarray.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/int_.hpp"
#include "pythonic/utils/numpy_conversion.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace details
  {
    template <class T, class U>
    long searchsorted(U const &a, T const &v, bool left)
    {
      if (left)
        return std::lower_bound(a.begin(), a.end(), v) - a.begin();
      else
        return std::upper_bound(a.begin(), a.end(), v) - a.begin();
    }

    bool issearchsortedleft(types::str const &side)
    {
      if (side[0] == "l")
        return true;
      else if (side[0] == "r")
        return false;
      else
        throw types::ValueError("'" + side + "' is an invalid value for keyword 'side'");
    }
  } // namespace details

  template <class T, class U>
  std::enable_if_t<!types::is_numexpr_arg<T>::value, long> searchsorted(U const &a, T const &v,
                                                                        types::str const &side)
  {
    bool left = details::issearchsortedleft(side);
    return details::searchsorted(a, v, left);
  }

  namespace
  {
    template <class E, class I0, class I1>
    void _search_sorted(E const &a, I0 ibegin, I0 iend, I1 obegin, bool left, utils::int_<1>)
    {
      for (; ibegin != iend; ++ibegin, ++obegin)
        *obegin = details::searchsorted(a, *ibegin, left);
    }

    template <class E, class I0, class I1, size_t N>
    void _search_sorted(E const &a, I0 ibegin, I0 iend, I1 obegin, bool left, utils::int_<N>)
    {
      for (; ibegin != iend; ++ibegin, ++obegin)
        _search_sorted(a, (*ibegin).begin(), (*ibegin).end(), (*obegin).begin(), left,
                       utils::int_<N - 1>());
    }
  } // namespace

  template <class E, class T>
  std::enable_if_t<types::is_numexpr_arg<E>::value,
                   types::ndarray<long, types::array_tuple<long, E::value>>>
  searchsorted(T const &a, E const &v, types::str const &side)
  {
    static_assert(T::value == 1, "Not Implemented : searchsorted for dimension != 1");
    bool left = details::issearchsortedleft(side);

    types::ndarray<long, types::array_tuple<long, E::value>> out(asarray(v)._shape, builtins::None);
    _search_sorted(a, v.begin(), v.end(), out.begin(), left, utils::int_<E::value>());
    return out;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
