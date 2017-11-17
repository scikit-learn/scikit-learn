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

float ATL_srefnrm2
(
   const int                  N,
   const float                * X,
   const int                  INCX
)
{
/*
 * Purpose
 * =======
 *
 * ATL_srefnrm2 returns the 2-norm of an n-vector x.
 *
 * Arguments
 * =========
 *
 * N       (input)                       const int
 *         On entry, N specifies the length of the vector x. N  must  be
 *         at least zero. Unchanged on exit.
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
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   register float             absxi, scale = ATL_sZERO,
                              ssq = ATL_sONE, temp, x0, x1, x2, x3, x4,
                              x5, x6, x7;
   float                      * StX = (float *)(X) + N * INCX;
   register int               i;
   int                        n = (int)(N), nu;
   const int                  incX2 = 2 * INCX, incX3 = 3 * INCX,
                              incX4 = 4 * INCX, incX5 = 5 * INCX,
                              incX6 = 6 * INCX, incX7 = 7 * INCX,
                              incX8 = 8 * INCX;
/* ..
 * .. Executable Statements ..
 *
 */
   if(    ( N  < 1 ) || ( INCX < 1 ) ) return( ATL_sZERO );
   else if( N == 1 )                   return( Msabs( *X ) );

   while( (X != StX)  && ( *X == ATL_sZERO ) ) { X += INCX; n--; }

   if( X == StX ) return( ATL_sZERO );

   if( ( nu = ( n >> 3 ) << 3 ) != 0 )
   {
      StX = (float *)X + nu * INCX;

      do
      {
         x0 = (*X);     x4 = X[incX4]; x1 = X[INCX ]; x5 = X[incX5];
         x2 = X[incX2]; x6 = X[incX6]; x3 = X[incX3]; x7 = X[incX7];

         absxi = Msabs( x0 );
         if( scale < absxi )
         {
            temp = scale / absxi; scale = absxi;
            ssq  = ATL_sONE + ssq * ( temp * temp );
         }
         else { temp = absxi / scale; ssq += temp * temp; }

         absxi = Msabs( x4 );
         if( scale < absxi )
         {
            temp = scale / absxi; scale = absxi;
            ssq  = ATL_sONE + ssq * ( temp * temp );
         }
         else { temp = absxi / scale; ssq += temp * temp; }

         absxi = Msabs( x1 );
         if( scale < absxi )
         {
            temp = scale / absxi; scale = absxi;
            ssq  = ATL_sONE + ssq * ( temp * temp );
         }
         else { temp = absxi / scale; ssq += temp * temp; }

         absxi = Msabs( x5 );
         if( scale < absxi )
         {
            temp = scale / absxi; scale = absxi;
            ssq  = ATL_sONE + ssq * ( temp * temp );
         }
         else { temp = absxi / scale; ssq += temp * temp; }

         absxi = Msabs( x2 );
         if( scale < absxi )
         {
            temp = scale / absxi; scale = absxi;
            ssq  = ATL_sONE + ssq * ( temp * temp );
         }
         else { temp = absxi / scale; ssq += temp * temp; }

         absxi = Msabs( x6 );
         if( scale < absxi )
         {
            temp = scale / absxi; scale = absxi;
            ssq  = ATL_sONE + ssq * ( temp * temp );
         }
         else { temp = absxi / scale; ssq += temp * temp; }

         absxi = Msabs( x3 );
         if( scale < absxi )
         {
            temp = scale / absxi; scale = absxi;
            ssq  = ATL_sONE + ssq * ( temp * temp );
         }
         else { temp = absxi / scale; ssq += temp * temp; }

         absxi = Msabs( x7 );
         if( scale < absxi )
         {
            temp = scale / absxi; scale = absxi;
            ssq  = ATL_sONE + ssq * ( temp * temp );
         }
         else { temp = absxi / scale; ssq += temp * temp; }

         X  += incX8;

      } while( X != StX );
   }

   for( i = n - nu; i != 0; i-- )
   {
      x0  = (*X);

      absxi = Msabs( x0 );
      if( scale < absxi )
      {
         temp = scale / absxi; scale = absxi;
         ssq  = ATL_sONE + ssq * ( temp * temp );
      }
      else { temp = absxi / scale; ssq += temp * temp; }

      X  += INCX;
   }
   return( scale * (float)(sqrt( (double)(ssq) )) );
/*
 * End of ATL_srefnrm2
 */
}
