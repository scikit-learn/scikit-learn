#include "blas.h"

int daxpy_(int *n, double *sa, double *sx, int *incx, double *sy,
           int *incy)
{
  long int i, m, ix, iy, nn, iincx, iincy;
  register double ssa;

  /* constant times a vector plus a vector.   
     uses unrolled loop for increments equal to one.   
     jack dongarra, linpack, 3/11/78.   
     modified 12/3/93, array(1) declarations changed to array(*) */

  /* Dereference inputs */
  nn = *n;
  ssa = *sa;
  iincx = *incx;
  iincy = *incy;

  if( nn > 0 && ssa != 0.0 )
  {
    if (iincx == 1 && iincy == 1) /* code for both increments equal to 1 */
    {
      m = nn-3;
      for (i = 0; i < m; i += 4)
      {
        sy[i] += ssa * sx[i];
        sy[i+1] += ssa * sx[i+1];
        sy[i+2] += ssa * sx[i+2];
        sy[i+3] += ssa * sx[i+3];
      }
      for ( ; i < nn; ++i) /* clean-up loop */
        sy[i] += ssa * sx[i];
    }
    else /* code for unequal increments or equal increments not equal to 1 */
    {
      ix = iincx >= 0 ? 0 : (1 - nn) * iincx;
      iy = iincy >= 0 ? 0 : (1 - nn) * iincy;
      for (i = 0; i < nn; i++)
      {
        sy[iy] += ssa * sx[ix];
        ix += iincx;
        iy += iincy;
      }
    }
  }

  return 0;
} /* daxpy_ */
