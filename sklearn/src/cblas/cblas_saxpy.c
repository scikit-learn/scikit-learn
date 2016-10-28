/* ---------------------------------------------------------------------
 *
 * -- Automatically Tuned Linear Algebra Software (ATLAS)
 *    (C) Copyright 2000 All Rights Reserved
 *
 * -- ATLAS routine -- Version 3.9.24 -- December 25, 2000
 *
 * Author         : Antoine P. Petitet
 * Originally developed at the University of Tennessee,
 * Innovative Computing Laboratory, Knoxville TN, 37996-1301, USA.
 *
 * ---------------------------------------------------------------------
 *
 * -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 * 1. Redistributions  of  source  code  must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce  the above copyright
 *    notice,  this list of conditions, and the  following disclaimer in
 *    the documentation and/or other materials provided with the distri-
 *    bution.
 * 3. The name of the University,  the ATLAS group,  or the names of its
 *    contributors  may not be used to endorse or promote products deri-
 *    ved from this software without specific written permission.
 *
 * -- Disclaimer:
 *
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,  INDIRECT, INCIDENTAL, SPE-
 * CIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO,  PROCUREMENT  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEO-
 * RY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  (IN-
 * CLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ---------------------------------------------------------------------
 */
/*
 * Include files
 */
#include "atlas_refmisc.h"

void cblas_saxpy
(
   const int                  N,
   const float                ALPHA,
   const float                * X,
   const int                  INCX,
   float                      * Y,
   const int                  INCY
)
{
/*
 * Purpose
 * =======
 *
 * ATL_srefaxpy performs the following operation:
 *
 *    y := y + alpha * x,
 *
 * where alpha is a scalar and x and y are two n-vectors.
 *
 * Arguments
 * =========
 *
 * N       (input)                       const int
 *         On entry, N specifies the length of the vector x. N  must  be
 *         at least zero. Unchanged on exit.
 *
 * ALPHA   (input)                       const float
 *         On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
 *         supplied as zero, then the entries of the incremented array X
 *         need not be set on input. Unchanged on exit.
 *
 * X       (input)                       const float *
 *         On entry,  X  points to the  first entry to be accessed of an
 *         incremented array of size equal to or greater than
 *            ( 1 + ( n - 1 ) * abs( INCX ) ) * sizeof(   float   ),
 *         that contains the vector x. Unchanged on exit.
 *
 * INCX    (input)                       const int
 *         On entry, INCX specifies the increment for the elements of X.
 *         INCX must not be zero. Unchanged on exit.
 *
 * Y       (input/output)                float *
 *         On entry,  Y  points to the  first entry to be accessed of an
 *         incremented array of size equal to or greater than
 *            ( 1 + ( n - 1 ) * abs( INCY ) ) * sizeof(   float   ),
 *         that contains the vector y.  On exit,  the entries of the in-
 *         cremented array  Y are updated with the scaled entries of the
 *         incremented array  X.
 *
 * INCY    (input)                       const int
 *         On entry, INCY specifies the increment for the elements of Y.
 *         INCY must not be zero. Unchanged on exit.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   register const float       alpha = ALPHA;
   register float             x0, x1, x2, x3, y0, y1, y2, y3;
   float                      * StX;
   register int               i;
   int                        nu;
   const int                  incX2 = 2 * INCX, incY2 = 2 * INCY,
                              incX3 = 3 * INCX, incY3 = 3 * INCY,
                              incX4 = 4 * INCX, incY4 = 4 * INCY;
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( N > 0 ) && ( alpha != ATL_sZERO ) )
   {
      if( ( nu = ( N >> 2 ) << 2 ) != 0 )
      {
         StX = (float *)X + nu * INCX;

         do
         {
            x0 = (*X);     y0 = (*Y);     x1 = X[INCX ]; y1 = Y[INCY ];
            x2 = X[incX2]; y2 = Y[incY2]; x3 = X[incX3]; y3 = Y[incY3];

            *Y       = y0 + alpha * x0; Y[INCY ] = y1 + alpha * x1;
            Y[incY2] = y2 + alpha * x2; Y[incY3] = y3 + alpha * x3;

            X  += incX4;
            Y  += incY4;

         } while( X != StX );
      }

      for( i = N - nu; i != 0; i-- )
      {
         x0  = (*X);
         y0  = (*Y);

         *Y  = y0 + alpha * x0;

         X  += INCX;
         Y  += INCY;
      }
   }
/*
 * End of ATL_srefaxpy
 */
}
