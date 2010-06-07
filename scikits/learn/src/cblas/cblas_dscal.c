/* ---------------------------------------------------------------------
 *
 * -- Automatically Tuned Linear Algebra Software (ATLAS)
 *    (C) Copyright 2000 All Rights Reserved
 *
 * -- ATLAS routine -- Version 3.2 -- December 25, 2000
 *
 * -- Suggestions,  comments,  bugs reports should be sent to the follo-
 *    wing e-mail address: atlas@cs.utk.edu
 *
 * Author         : Antoine P. Petitet
 * University of Tennessee - Innovative Computing Laboratory
 * Knoxville TN, 37996-1301, USA.
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

void cblas_dscal
(
   const int                  N,
   const double               ALPHA,
   double                     * X,
   const int                  INCX
)
{
/*
 * Purpose
 * =======
 *
 * ATL_drefscal performs the following operation:
 *
 *    x := alpha * x,
 *
 * where alpha is a scalar and x is an n-vector.
 *
 * Arguments
 * =========
 *
 * N       (input)                       const int
 *         On entry, N specifies the length of the vector x. N  must  be
 *         at least zero. Unchanged on exit.
 *
 * ALPHA   (input)                       const double
 *         On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
 *         supplied as zero, then the entries of the incremented array X
 *         need not be set on input. Unchanged on exit.
 *
 * X       (input/output)                double *
 *         On entry,  X  points to the  first entry to be accessed of an
 *         incremented array of size equal to or greater than
 *            ( 1 + ( n - 1 ) * abs( INCX ) ) * sizeof(   double  ),
 *         that contains the vector x.  On exit,  the entries of the in-
 *         cremented array X are mutiplied by alpha.
 *
 * INCX    (input)                       const int
 *         On entry, INCX specifies the increment for the elements of X.
 *         INCX must not be zero. Unchanged on exit.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
/* ..
 * .. Executable Statements ..
 *
 */
   register double            x0, x1, x2, x3, x4, x5, x6, x7;
   register const double      alpha = ALPHA;
   double                     * StX;
   register int               i;
   int                        nu;
   const int                  incX2 = 2 * INCX, incX3 = 3 * INCX,
                              incX4 = 4 * INCX, incX5 = 5 * INCX,
                              incX6 = 6 * INCX, incX7 = 7 * INCX,
                              incX8 = 8 * INCX;
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( N > 0 ) && ( alpha != ATL_dONE ) )
   {
      if( alpha == ATL_dZERO )
      {
         if( ( nu = ( N >> 3 ) << 3 ) != 0 )
         {
            StX = (double *)X + nu * INCX;

            do
            {
               (*X)     = ATL_dZERO; X[incX4] = ATL_dZERO;
               X[INCX ] = ATL_dZERO; X[incX5] = ATL_dZERO;
               X[incX2] = ATL_dZERO; X[incX6] = ATL_dZERO;
               X[incX3] = ATL_dZERO; X[incX7] = ATL_dZERO;

               X  += incX8;

            } while( X != StX );
         }

         for( i = N - nu; i != 0; i-- )
         {
            *X = ATL_dZERO;
            X += INCX;
         }
      }
      else
      {
         if( ( nu = ( N >> 3 ) << 3 ) != 0 )
         {
            StX = (double *)X + nu * INCX;

            do
            {
               x0 = (*X);     x4 = X[incX4];
               x1 = X[INCX ]; x5 = X[incX5];
               x2 = X[incX2]; x6 = X[incX6];
               x3 = X[incX3]; x7 = X[incX7];

               x0 *= alpha;   x4 *= alpha;
               x1 *= alpha;   x5 *= alpha;
               x2 *= alpha;   x6 *= alpha;
               x3 *= alpha;   x7 *= alpha;

               (*X)     = x0; X[incX4] = x4;
               X[INCX ] = x1; X[incX5] = x5;
               X[incX2] = x2; X[incX6] = x6;
               X[incX3] = x3; X[incX7] = x7;

               X  += incX8;

            } while( X != StX );
         }

         for( i = N - nu; i != 0; i-- )
         {
            x0  = (*X);
            x0 *= alpha;
            *X  = x0;
            X  += INCX;
         }
      }
   }
/*
 * End of ATL_drefscal
 */
}
