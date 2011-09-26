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
#include "atlas_reflvl2.h"
#include "atlas_reflevel2.h"

void ATL_sreftrsv
(
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const enum ATLAS_DIAG      DIAG,
   const int                  N,
   const float                * A,
   const int                  LDA,
   float                      * X,
   const int                  INCX
)
{
/*
 * Purpose
 * =======
 *
 * ATL_sreftrsv solves one of the systems of equations
 *
 *    A * x = b,   or   A' * x = b,
 *
 * where b and x are n-element vectors and  A is an n by n unit, or non-
 * unit, upper or lower triangular matrix.
 *
 * No test for  singularity  or  near-singularity  is included  in  this
 * routine. Such tests must be performed before calling this routine.
 *
 * Arguments
 * =========
 *
 * UPLO    (input)                       const enum ATLAS_UPLO
 *         On entry, UPLO  specifies whether  the  matrix is an upper or
 *         lower triangular matrix as follows:
 *
 *             UPLO = AtlasUpper   A is an upper triangular matrix.
 *
 *             UPLO = AtlasLower   A is a lower triangular matrix.
 *
 *         Unchanged on exit.
 *
 * TRANS   (input)                       const enum ATLAS_TRANS
 *         On entry,  TRANS specifies the equations to be solved as fol-
 *         lows:
 *
 *            TRANS = AtlasNoTrans     A  * x = b,
 *
 *            TRANS = AtlasConj        A  * x = b,
 *
 *            TRANS = AtlasTrans       A' * x = b,
 *
 *            TRANS = AtlasTrans       A' * x = b.
 *
 *         Unchanged on exit.
 *
 * DIAG    (input)                       const enum ATLAS_DIAG
 *         On entry, DIAG specifies whether or not A is unit triangu-
 *         lar as follows:
 *
 *            DIAG = AtlasUnit       A is assumed to be unit triangular,
 *
 *            DIAG = AtlasNonUnit    A is not assumed to be unit trian-
 *                                   gular.
 *
 *         Unchanged on exit.
 *
 * N       (input)                       const int
 *         On entry, N specifies the order of the matrix A. N must be at
 *         least zero. Unchanged on exit.
 *
 * A       (input)                       const float *
 *         On entry,  A  points  to an array of size equal to or greater
 *         than   LDA * n * sizeof(   float   ).   Before   entry   with
 *         UPLO = AtlasUpper, the leading  n by n  upper triangular part
 *         of the array  A must contain the upper triangular  matrix and
 *         the strictly lower triangular part of  A  is  not referenced.
 *         Before entry with UPLO = AtlasLower, the leading n by n lower
 *         triangular  part of the array A must contain the lower trian-
 *         gular matrix  and the strictly upper triangular part  of A is
 *         not referenced. Unchanged on exit.
 *
 *         Note that when  DIAG = AtlasUnit,  the diagonal elements of A
 *         are not referenced  either,  but are assumed to be unity.
 *
 * LDA     (input)                       const int
 *         On entry, LDA  specifies the leading dimension of A as decla-
 *         red  in  the  calling  (sub) program.  LDA  must be  at least
 *         MAX( 1, n ). Unchanged on exit.
 *
 * X       (input/output)                float *
 *         On entry,  X  points to the  first entry to be accessed of an
 *         incremented array of size equal to or greater than
 *            ( 1 + ( n - 1 ) * abs( INCX ) ) * sizeof(   float   ),
 *         that contains the vector x. Before entry, the incremented ar-
 *         ray X must contain the n element right-hand side vector b. On
 *         exit, X is overwritten with the solution vector x.
 *
 * INCX    (input)                       const int
 *         On entry, INCX specifies the increment for the elements of X.
 *         INCX must not be zero. Unchanged on exit.
 *
 * ---------------------------------------------------------------------
 */
/* ..
 * .. Executable Statements ..
 *
 */
   if( N == 0 ) return;

   if( UPLO == AtlasUpper )
   {
      if( ( TRANS == AtlasNoTrans ) || ( TRANS == AtlasConj ) )
      {
         if( DIAG == AtlasNonUnit )
         {
            ATL_sreftrsvUNN( N,    A, LDA, X, INCX );
         }
         else
         {
            ATL_sreftrsvUNU( N,    A, LDA, X, INCX );
         }
      }
      else
      {
         if( DIAG == AtlasNonUnit )
         {
            ATL_sreftrsvUTN( N,    A, LDA, X, INCX );
         }
         else
         {
            ATL_sreftrsvUTU( N,    A, LDA, X, INCX );
         }
      }
   }
   else
   {
      if( ( TRANS == AtlasNoTrans ) || ( TRANS == AtlasConj ) )
      {
         if( DIAG == AtlasNonUnit )
         {
            ATL_sreftrsvLNN( N,    A, LDA, X, INCX );
         }
         else
         {
            ATL_sreftrsvLNU( N,    A, LDA, X, INCX );
         }
      }
      else
      {
         if( DIAG == AtlasNonUnit )
         {
            ATL_sreftrsvLTN( N,    A, LDA, X, INCX );
         }
         else
         {
            ATL_sreftrsvLTU( N,    A, LDA, X, INCX );
         }
      }
   }
/*
 * End of ATL_sreftrsv
 */
}
