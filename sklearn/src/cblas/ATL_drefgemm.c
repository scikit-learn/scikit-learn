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
#include "atlas_reflvl3.h"
#include "atlas_reflevel3.h"

void ATL_drefgemm
(
   const enum ATLAS_TRANS     TRANSA,
   const enum ATLAS_TRANS     TRANSB,
   const int                  M,
   const int                  N,
   const int                  K,
   const double               ALPHA,
   const double               * A,
   const int                  LDA,
   const double               * B,
   const int                  LDB,
   const double               BETA,
   double                     * C,
   const int                  LDC
)
{
/*
 * Purpose
 * =======
 *
 * ATL_drefgemm  performs one of the matrix-matrix operations
 *
 *    C := alpha * op( A ) * op( B ) + beta * C,
 *
 * where op( X ) is one of
 *
 *    op( X ) = X   or   op( X ) = X'.
 *
 * Alpha and beta are scalars, and A, B and C are matrices, with op( A )
 * an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
 *
 * Arguments
 * =========
 *
 * TRANSA  (input)                       const enum ATLAS_TRANS
 *         On entry, TRANSA  specifies the form of op( A ) to be used in
 *         the matrix multiplication as follows:
 *
 *            TRANSA = AtlasNoTrans    op( A ) = A,
 *
 *            TRANSA = AtlasTrans      op( A ) = A',
 *
 *            TRANSA = AtlasConjTrans  op( A ) = A'.
 *
 *         Unchanged on exit.
 *
 * TRANSB  (input)                       const enum ATLAS_TRANS
 *         On entry, TRANSB  specifies the form of op( A ) to be used in
 *         the matrix multiplication as follows:
 *
 *            TRANSB = AtlasNoTrans    op( B ) = B,
 *
 *            TRANSB = AtlasTrans      op( B ) = B',
 *
 *            TRANSB = AtlasConjTrans  op( B ) = B'.
 *
 *         Unchanged on exit.
 *
 * M       (input)                       const int
 *         On entry,  M  specifies  the  number  of rows  of the  matrix
 *         op( A )  and  of the  matrix  C.  M  must  be at least  zero.
 *         Unchanged on exit.
 *
 * N       (input)                       const int
 *         On entry,  N  specifies  the number  of columns of the matrix
 *         op( B )  and the number of columns of the matrix C. N must be
 *         at least zero. Unchanged on exit.
 *
 * K       (input)                       const int
 *         On entry,  K  specifies  the  number of columns of the matrix
 *         op( A ) and the number of rows  of the matrix op( B ). K must
 *         be at least  zero. Unchanged on exit.
 *
 * ALPHA   (input)                       const double
 *         On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
 *         supplied  as  zero  then the elements of the matrices A and B
 *         need not be set on input. Unchanged on exit.
 *
 * A       (input)                       const double *
 *         On entry,  A  points  to an array of size equal to or greater
 *         than   LDA * ka * sizeof(   double  ),   where  ka  is k when
 *         TRANSA = AtlasNoTrans, and is m otherwise. Before  entry with
 *         TRANSA = AtlasNoTrans, the leading m by k part of the array A
 *         must contain the matrix  A, otherwise the leading k by m part
 *         of the array A must contain the matrix A. Unchanged on exit.
 *
 * LDA     (input)                       const int
 *         On entry, LDA  specifies the leading dimension of A as decla-
 *         red  in  the  calling  (sub) program.  LDA  must be  at least
 *         MAX( 1, m ) when TRANS = AtlasNotrans, and MAX( 1, k ) other-
 *         wise. Unchanged on exit.
 *
 * B       (input)                       const double *
 *         On entry,  B  points  to an array of size equal to or greater
 *         than   LDB * kb * sizeof(   double  ),   where  kb  is n when
 *         TRANSB = AtlasNoTrans, and is k otherwise. Before  entry with
 *         TRANSB = AtlasNoTrans, the leading k by n part of the array B
 *         must contain the matrix  B, otherwise the leading n by k part
 *         of the array B must contain the matrix B. Unchanged on exit.
 *
 * LDB     (input)                       const int
 *         On entry, LDB  specifies the leading dimension of A as decla-
 *         red  in  the  calling  (sub) program.  LDB  must be  at least
 *         MAX( 1, k )  when  TRANS = AtlasNotrans or TRANS = AtlasConj,
 *         and MAX( 1, n ) otherwise. Unchanged on exit.
 *
 * BETA    (input)                       const double
 *         On entry,  BETA  specifies the scalar  beta.   When  BETA  is
 *         supplied  as  zero  then  the  elements of the matrix C  need
 *         not be set on input. Unchanged on exit.
 *
 * C       (input/output)                double *
 *         On entry,  C  points  to an array of size equal to or greater
 *         than   LDC * n * sizeof(   double  ). Before  entry, the lea-
 *         ding  m by n  part of the array C must contain the matrix  C,
 *         except when beta is zero, in which case C need not be  set on
 *         entry. On exit, the array C is overwritten by the  m by n ma-
 *         trix ( alpha*op( A )*op( B ) + beta*C ).
 *
 * LDC     (input)                       const int
 *         On entry, LDC  specifies the leading dimension of A as decla-
 *         red  in  the  calling  (sub) program.  LDC  must be  at least
 *         MAX( 1, m ). Unchanged on exit.
 *
 * ---------------------------------------------------------------------
 */
   if( ( M == 0 ) || ( N == 0 ) ||
       ( ( ( ALPHA == ATL_dZERO ) || ( K == 0 ) ) &&
         ( BETA == ATL_dONE ) ) ) return;

   if( ALPHA == ATL_dZERO )
   { Mdgescal( M, N, BETA, C, LDC ); return; }

   if( TRANSB == AtlasNoTrans )
   {
      if( TRANSA == AtlasNoTrans )
      {
         ATL_drefgemmNN( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC );
      }
      else
      {
         ATL_drefgemmTN( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC );
      }
   }
   else
   {
      if( TRANSA == AtlasNoTrans )
      {
         ATL_drefgemmNT( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC );
      }
      else
      {
         ATL_drefgemmTT( M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC );
      }
   }
/*
 * End of ATL_drefgemm
 */
}
