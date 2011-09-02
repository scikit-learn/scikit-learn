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
#include "atlas_reflevel1.h"

void ATL_drefrot
(
   const int                  N,
   double                     * X,
   const int                  INCX,
   double                     * Y,
   const int                  INCY,
   const double               C,
   const double               S
)
{
/*
 * Purpose
 * =======
 *
 * ATL_drefrot  applies a  plane  rotation to the two n-vectors x and y.
 * This routine computes:
 *
 *    [ x_i ]   [ c  s ] [ x_i ]
 *    [ y_i ] = [ -s c ] [ y_i ]    for all i = 1 .. n.
 *
 * If n <= 0 or if c = 1 and s = 0, this subroutine returns immediately.
 *
 * Arguments
 * =========
 *
 * N       (input)                       const int
 *         On entry, N specifies the length of the vector x. N  must  be
 *         at least zero. Unchanged on exit.
 *
 * X       (input/output)                double *
 *         On entry,  X  points to the  first entry to be accessed of an
 *         incremented array of size equal to or greater than
 *            ( 1 + ( n - 1 ) * abs( INCX ) ) * sizeof(   double  ),
 *         that contains the vector x.  On exit,  the entries of the in-
 *         cremented array  X are rotated with the entries of the incre-
 *         mented array  Y.
 *
 * INCX    (input)                       const int
 *         On entry, INCX specifies the increment for the elements of X.
 *         INCX must not be zero. Unchanged on exit.
 *
 * Y       (input/output)                double *
 *         On entry,  Y  points to the  first entry to be accessed of an
 *         incremented array of size equal to or greater than
 *            ( 1 + ( n - 1 ) * abs( INCY ) ) * sizeof(   double  ),
 *         that contains the vector y.  On exit,  the entries of the in-
 *         cremented array  Y are rotated with the entries of the incre-
 *         mented array  X.
 *
 * INCY    (input)                       const int
 *         On entry, INCY specifies the increment for the elements of Y.
 *         INCY must not be zero. Unchanged on exit.
 *
 * C       (input)                       const double
 *         On entry, C specifies the scalar c definiting the plane rota-
 *         tion. Unchanged on exit.
 *
 * S       (input)                       const double
 *         On entry, S specifies the scalar s definiting the plane rota-
 *         tion. Unchanged on exit.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   register double            x0, x1, x2, x3, y0, y1, y2, y3;
   register const double      co = C, si = S;
   double                     * StX;
   register int               i;
   int                        nu;
   const int                  incX2 = 2 * INCX, incY2 = 2 * INCY,
                              incX3 = 3 * INCX, incY3 = 3 * INCY,
                              incX4 = 4 * INCX, incY4 = 4 * INCY;
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( N > 0 ) && !Mdone( co, si ) )
   {
      if( ( nu = ( N >> 2 ) << 2 ) != 0 )
      {
         StX = (double *)X + nu * INCX;

         do
         {
            x0 = (*X);     y0 = (*Y);
            x1 = X[INCX ]; y1 = Y[INCY ];
            x2 = X[incX2]; y2 = Y[incY2];
            x3 = X[incX3]; y3 = Y[incY3];

            *X       = co * x0 + si * y0; *Y       = co * y0 - si * x0;
            X[INCX ] = co * x1 + si * y1; Y[INCY ] = co * y1 - si * x1;
            X[incX2] = co * x2 + si * y2; Y[incY2] = co * y2 - si * x2;
            X[incX3] = co * x3 + si * y3; Y[incY3] = co * y3 - si * x3;

            X  += incX4;
            Y  += incY4;

         } while( X != StX );
      }

      for( i = N - nu; i != 0; i-- )
      {
         x0  = (*X);
         y0  = (*Y);

         *X  = co * x0 + si * y0;
         *Y  = co * y0 - si * x0;

         X  += INCX;
         Y  += INCY;
      }
   }
/*
 * End of ATL_drefrot
 */
}
