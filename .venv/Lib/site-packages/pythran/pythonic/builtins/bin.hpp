#ifndef PYTHONIC_BUILTIN_BIN_HPP
#define PYTHONIC_BUILTIN_BIN_HPP

#include "pythonic/include/builtins/bin.hpp"

#include "pythonic/types/str.hpp"
#include "pythonic/utils/functor.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

PYTHONIC_NS_BEGIN

namespace builtins
{
  template <class T>
  std::enable_if_t<std::is_scalar<T>::value, types::str> bin(T const &v)
  {
    using UT = typename std::make_unsigned<T>::type;
    if (v < T{0})
      if (v == std::numeric_limits<T>::min()) {
        // In this special case, calling -v would overflow so
        // a special case is needed.
        types::str res;
        auto &backend = res.chars();
        backend.resize(8 * sizeof(T) + 3);
        auto it = backend.begin();
        *it++ = '-';
        *it++ = '0';
        *it++ = 'b';
        *it++ = '1';
        std::fill(it, backend.end(), '0');
        return res;
      } else
        return "-" + bin(-v);
    else if (v == T{0})
      return "0b0";
    else {
      // Due to rounding errors, we cannot use std::log2(v)
      // to accuratly find length.
      size_t len = (8 * sizeof(UT)) - 1;
      UT i{UT{1} << len};
      while (!(i & v)) {
        i >>= 1;
        len--;
      }
      types::str res;
      res.reserve(2 + len);
      auto &backend = res.chars();
      backend.append("0b");
      for (; i; i >>= 1)
        if (v & i)
          backend.append("1");
        else
          backend.append("0");
      return res;
    }
  }
} // namespace builtins
PYTHONIC_NS_END

#endif
