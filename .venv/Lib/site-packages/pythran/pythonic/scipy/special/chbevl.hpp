
/* Adapated from Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1985, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

#ifndef PYTHONIC_SCIPY_SPECIAL_CHBEVL_HPP
#define PYTHONIC_SCIPY_SPECIAL_CHBEVL_HPP

PYTHONIC_NS_BEGIN

namespace scipy
{
  namespace special
  {
    template <size_t N>
    double chbevl(double x, const double (&array)[N])
    {
      const double *p = &array[0];
      double b0 = *p++;
      double b1 = 0.0;
      double b2;
      size_t i = N - 1;

      do {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + *p++;
      } while (--i);

      return (0.5 * (b0 - b2));
    }

  } // namespace special
} // namespace scipy
PYTHONIC_NS_END

#endif
