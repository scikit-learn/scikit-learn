#ifndef PYTHONIC_NUMPY_BASEREPR_HPP
#define PYTHONIC_NUMPY_BASEREPR_HPP

#include "pythonic/include/numpy/base_repr.hpp"

#include "pythonic/types/ndarray.hpp"
#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  types::str base_repr(long number, long base, long padding)
  {
    types::str res;

    // check that the base if valid
    if (base < 2 || base > 16) {
      return res;
    }

    int const ndigits = (number == 0 ? 1 : std::ceil(std::log(std::labs(number)) / std::log(base)));
    int const effective_padding = padding - ((number == 0) && (padding > 0) ? 1 : 0);

    res.resize(ndigits + effective_padding + (number < 0 ? 1 : 0));

    // Apply negative sign
    auto it = res.chars().begin();
    if (number < 0)
      *it++ = '-';

    // Apply padding
    std::fill(it, std::next(it, effective_padding), '0');

    auto rit = res.chars().rbegin();
    long quotient = std::labs(number);

    do {
      const long tmp = quotient / base;
      *rit++ = "0123456789ABCDEF"[quotient - (tmp * base)];
      quotient = tmp;
    } while (quotient);

    return res;
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
