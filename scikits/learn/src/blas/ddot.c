#include "blas.h"

double ddot_(int *n, double *sx, int *incx, double *sy, int *incy)
{
  long int i, m, nn, iincx, iincy;
  double stemp;
  long int ix, iy;

  /* forms the dot product of two vectors.   
     uses unrolled loops for increments equal to one.   
     jack dongarra, linpack, 3/11/78.   
     modified 12/3/93, array(1) declarations changed to array(*) */

  /* Dereference inputs */
  nn = *n;
  iincx = *incx;
  iincy = *incy;

  stemp = 0.0;
  if (nn > 0)
  {
    if (iincx == 1 && iincy == 1) /* code for both increments equal to 1 */
    {
      m = nn-4;
      for (i = 0; i < m; i += 5)
        stemp += sx[i] * sy[i] + sx[i+1] * sy[i+1] + sx[i+2] * sy[i+2] +
                 sx[i+3] * sy[i+3] + sx[i+4] * sy[i+4];

      for ( ; i < nn; i++)        /* clean-up loop */
        stemp += sx[i] * sy[i];
    }
    else /* code for unequal increments or equal increments not equal to 1 */
    {
      ix = 0;
      iy = 0;
      if (iincx < 0)
        ix = (1 - nn) * iincx;
      if (iincy < 0)
        iy = (1 - nn) * iincy;
      for (i = 0; i < nn; i++)
      {
        stemp += sx[ix] * sy[iy];
        ix += iincx;
        iy += iincy;
      }
    }
  }

  return stemp;
} /* ddot_ */
