/*
 *             Automatically Tuned Linear Algebra Software v3.9.25
 *                    (C) Copyright 1999 R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */
#ifndef ATLAS_LEVEL2_H
#define ATLAS_LEVEL2_H
#include "atlas_refalias2.h"

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void ATL_sgemv(const enum ATLAS_TRANS TransA, const int M, const int N,
               const float alpha, const float *A, const int lda,
               const float *X, const int incX, const float beta,
               float *Y, const int incY);
void ATL_sgemvT_L2
   (ATL_CINT M, ATL_CINT N, const float alpha, const float *A, ATL_CINT lda,
    const float *X, ATL_CINT incX, const float beta, float *Y, ATL_CINT incY);
void ATL_sgemvT_L1
   (ATL_CINT M, ATL_CINT N, const float alpha, const float *A, ATL_CINT lda,
    const float *X, ATL_CINT incX, const float beta, float *Y, ATL_CINT incY);
void ATL_sgemvT
   (ATL_CINT M, ATL_CINT N, const float alpha, const float *A, ATL_CINT lda,
    const float *X, ATL_CINT incX, const float beta, float *Y, ATL_CINT incY);
void ATL_sgbmv(const enum ATLAS_TRANS TransA, const int M, const int N,
               const int KL, const int KU, const float alpha,
               const float *A, const int lda, const float *X,
               const int incX, const float beta, float *Y, const int incY);
void ATL_strmv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N,
               const float *A, const int lda, float *X, const int incX);
void ATL_stbmv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N, const int K,
               const float *A, const int lda, float *X, const int incX);
void ATL_stpmv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N, const float *Ap,
               float *X, const int incX);
void ATL_strsv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N,
               const float *A, const int lda, float *X, const int incX);
void ATL_stbsv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N, const int K,
               const float *A, const int lda, float *X, const int incX);
void ATL_stpsv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N,
               const float *Ap, float *X, const int incX);

void ATL_dgemv(const enum ATLAS_TRANS TransA, const int M, const int N,
               const double alpha, const double *A, const int lda,
               const double *X, const int incX, const double beta,
               double *Y, const int incY);
void ATL_dgemvT_L2
   (ATL_CINT M, ATL_CINT N, const double alpha, const double *A, ATL_CINT lda,
    const double *X, ATL_CINT incX, const double beta, double *Y, ATL_CINT incY);
void ATL_dgemvT_L1
   (ATL_CINT M, ATL_CINT N, const double alpha, const double *A, ATL_CINT lda,
    const double *X, ATL_CINT incX, const double beta, double *Y, ATL_CINT incY);
void ATL_dgemvT
   (ATL_CINT M, ATL_CINT N, const double alpha, const double *A, ATL_CINT lda,
    const double *X, ATL_CINT incX, const double beta, double *Y, ATL_CINT incY);
void ATL_dgbmv(const enum ATLAS_TRANS TransA, const int M, const int N,
               const int KL, const int KU, const double alpha,
               const double *A, const int lda, const double *X,
               const int incX, const double beta, double *Y, const int incY);
void ATL_dtrmv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N,
               const double *A, const int lda, double *X, const int incX);
void ATL_dtbmv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N, const int K,
               const double *A, const int lda, double *X, const int incX);
void ATL_dtpmv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N, const double *Ap,
               double *X, const int incX);
void ATL_dtrsv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N,
               const double *A, const int lda, double *X, const int incX);
void ATL_dtbsv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N, const int K,
               const double *A, const int lda, double *X, const int incX);
void ATL_dtpsv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N,
               const double *Ap, double *X, const int incX);

void ATL_cgemv(const enum ATLAS_TRANS TransA, const int M, const int N,
               const float *alpha, const float *A, const int lda,
               const float *X, const int incX, const float *beta,
               float *Y, const int incY);
void ATL_cgemvT_L2
   (ATL_CINT M, ATL_CINT N, const float *alpha, const float *A, ATL_CINT lda,
    const float *X, ATL_CINT incX, const float *beta, float *Y, ATL_CINT incY);
void ATL_cgemvT_L1
   (ATL_CINT M, ATL_CINT N, const float *alpha, const float *A, ATL_CINT lda,
    const float *X, ATL_CINT incX, const float *beta, float *Y, ATL_CINT incY);
void ATL_cgemvT
   (ATL_CINT M, ATL_CINT N, const float *alpha, const float *A, ATL_CINT lda,
    const float *X, ATL_CINT incX, const float *beta, float *Y, ATL_CINT incY);
void ATL_cgbmv(const enum ATLAS_TRANS TransA, const int M, const int N,
               const int KL, const int KU, const float *alpha,
               const float *A, const int lda, const float *X,
               const int incX, const float *beta, float *Y, const int incY);
void ATL_ctrmv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N,
               const float *A, const int lda, float *X, const int incX);
void ATL_ctbmv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N, const int K,
               const float *A, const int lda, float *X, const int incX);
void ATL_ctpmv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N, const float *Ap,
               float *X, const int incX);
void ATL_ctrsv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N,
               const float *A, const int lda, float *X, const int incX);
void ATL_ctbsv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N, const int K,
               const float *A, const int lda, float *X, const int incX);
void ATL_ctpsv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N,
               const float *Ap, float *X, const int incX);

void ATL_zgemv(const enum ATLAS_TRANS TransA, const int M, const int N,
               const double *alpha, const double *A, const int lda,
               const double *X, const int incX, const double *beta,
               double *Y, const int incY);
void ATL_zgemvT_L2
   (ATL_CINT M, ATL_CINT N, const double *alpha, const double *A, ATL_CINT lda,
    const double *X, ATL_CINT incX, const double *beta, double *Y, ATL_CINT incY);
void ATL_zgemvT_L1
   (ATL_CINT M, ATL_CINT N, const double *alpha, const double *A, ATL_CINT lda,
    const double *X, ATL_CINT incX, const double *beta, double *Y, ATL_CINT incY);
void ATL_zgemvT
   (ATL_CINT M, ATL_CINT N, const double *alpha, const double *A, ATL_CINT lda,
    const double *X, ATL_CINT incX, const double *beta, double *Y, ATL_CINT incY);
void ATL_zgbmv(const enum ATLAS_TRANS TransA, const int M, const int N,
               const int KL, const int KU, const double *alpha,
               const double *A, const int lda, const double *X,
               const int incX, const double *beta, double *Y, const int incY);
void ATL_ztrmv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N,
               const double *A, const int lda, double *X, const int incX);
void ATL_ztbmv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N, const int K,
               const double *A, const int lda, double *X, const int incX);
void ATL_ztpmv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N, const double *Ap,
               double *X, const int incX);
void ATL_ztrsv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N,
               const double *A, const int lda, double *X, const int incX);
void ATL_ztbsv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N, const int K,
               const double *A, const int lda, double *X, const int incX);
void ATL_ztpsv(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TransA,
               const enum ATLAS_DIAG Diag, const int N,
               const double *Ap, double *X, const int incX);


/*
 * Routines with S and D prefixes only
 */
void ATL_ssymv(const enum ATLAS_UPLO Uplo, const int N,
               const float alpha, const float *A, const int lda,
               const float *X, const int incX, const float beta,
               float *Y, const int incY);
void ATL_ssbmv(const enum ATLAS_UPLO Uplo, const int N, const int K,
               const float alpha, const float *A, const int lda,
               const float *X, const int incX, const float beta,
               float *Y, const int incY);
void ATL_sspmv(const enum ATLAS_UPLO Uplo, const int N, const float alpha,
               const float *Ap, const float *X, const int incX,
               const float beta, float *Y, const int incY);
void ATL_sger(const int M, const int N, const float alpha,
              const float *X, const int incX, const float *Y, const int incY,
              float *A, const int lda);
void ATL_sger2(const int M, const int N, const float alpha,
               const float *X, const int incX, const float *Y, const int incY,
               const float beta,
               const float *W, const int incW, const float *Z, const int incZ,
               float *A, const int lda);
void ATL_ssyr(const enum ATLAS_UPLO Uplo, const int N, const float alpha,
              const float *X, const int incX, float *A, const int lda);
void ATL_sspr(const enum ATLAS_UPLO Uplo, const int N, const float alpha,
              const float *X, const int incX, float *Ap);
void ATL_ssyr2(const enum ATLAS_UPLO Uplo, const int N, const float alpha,
               const float *X, const int incX, const float *Y, const int incY,
               float *A, const int lda);
void ATL_sspr2(const enum ATLAS_UPLO Uplo, const int N, const float alpha,
               const float *X, const int incX, const float *Y, const int incY,
               float *A);

void ATL_dsymv(const enum ATLAS_UPLO Uplo, const int N,
               const double alpha, const double *A, const int lda,
               const double *X, const int incX, const double beta,
               double *Y, const int incY);
void ATL_dsbmv(const enum ATLAS_UPLO Uplo, const int N, const int K,
               const double alpha, const double *A, const int lda,
               const double *X, const int incX, const double beta,
               double *Y, const int incY);
void ATL_dspmv(const enum ATLAS_UPLO Uplo, const int N, const double alpha,
               const double *Ap, const double *X, const int incX,
               const double beta, double *Y, const int incY);
void ATL_dger(const int M, const int N, const double alpha,
              const double *X, const int incX, const double *Y, const int incY,
              double *A, const int lda);
void ATL_dger2(const int M, const int N, const double alpha,
               const double *X, const int incX, const double *Y, const int incY,
               const double beta,
               const double *W, const int incW, const double *Z, const int incZ,
               double *A, const int lda);
void ATL_dsyr(const enum ATLAS_UPLO Uplo, const int N, const double alpha,
              const double *X, const int incX, double *A, const int lda);
void ATL_dspr(const enum ATLAS_UPLO Uplo, const int N, const double alpha,
              const double *X, const int incX, double *Ap);
void ATL_dsyr2(const enum ATLAS_UPLO Uplo, const int N, const double alpha,
               const double *X, const int incX, const double *Y, const int incY,
               double *A, const int lda);
void ATL_dspr2(const enum ATLAS_UPLO Uplo, const int N, const double alpha,
               const double *X, const int incX, const double *Y, const int incY,
               double *A);


/*
 * Routines with C and Z prefixes only
 */
void ATL_chemv(const enum ATLAS_UPLO Uplo, const int N,
               const float *alpha, const float *A, const int lda,
               const float *X, const int incX, const float *beta,
               float *Y, const int incY);
void ATL_chbmv(const enum ATLAS_UPLO Uplo, const int N, const int K,
               const float *alpha, const float *A, const int lda,
               const float *X, const int incX, const float *beta,
               float *Y, const int incY);
void ATL_chpmv(const enum ATLAS_UPLO Uplo, const int N,
               const float *alpha, const float *Ap,
               const float *X, const int incX, const float *beta,
               float *Y, const int incY);
void ATL_cgeru(const int M, const int N, const float *alpha,
               const float *X, const int incX, const float *Y, const int incY,
               float *A, const int lda);
void ATL_cgerc(const int M, const int N, const float *alpha,
               const float *X, const int incX, const float *Y, const int incY,
               float *A, const int lda);
void ATL_cger2u(const int M, const int N,
                const float *alpha, const float *X, const int incX,
                const float *Y, const int incY,
                const float *beta, const float *W, const int incW,
                const float *Z, const int incZ,
                float *A, const int lda);
void ATL_cger2c(const int M, const int N,
                const float *alpha, const float *X, const int incX,
                const float *Y, const int incY,
                const float *beta, const float *W, const int incW,
                const float *Z, const int incZ,
                float *A, const int lda);
void ATL_cher(const enum ATLAS_UPLO Uplo, const int N, const float alpha,
              const float *X, const int incX, float *A, const int lda);
void ATL_chpr(const enum ATLAS_UPLO Uplo, const int N, const float alpha,
                   const float *X, const int incX, float *A);
void ATL_cher2(const enum ATLAS_UPLO Uplo, const int N,
               const float *alpha, const float *X, const int incX,
               const float *Y, const int incY, float *A, const int lda);
void ATL_chpr2(const enum ATLAS_UPLO Uplo, const int N,
               const float *alpha, const float *X, const int incX,
               const float *Y, const int incY, float *Ap);

void ATL_zhemv(const enum ATLAS_UPLO Uplo, const int N,
               const double *alpha, const double *A, const int lda,
               const double *X, const int incX, const double *beta,
               double *Y, const int incY);
void ATL_zhbmv(const enum ATLAS_UPLO Uplo, const int N, const int K,
               const double *alpha, const double *A, const int lda,
               const double *X, const int incX, const double *beta,
               double *Y, const int incY);
void ATL_zhpmv(const enum ATLAS_UPLO Uplo, const int N,
               const double *alpha, const double *Ap,
               const double *X, const int incX, const double *beta,
               double *Y, const int incY);
void ATL_zgeru(const int M, const int N, const double *alpha,
               const double *X, const int incX, const double *Y, const int incY,
               double *A, const int lda);
void ATL_zgerc(const int M, const int N, const double *alpha,
               const double *X, const int incX, const double *Y, const int incY,
               double *A, const int lda);
void ATL_zger2u(const int M, const int N,
                const double *alpha, const double *X, const int incX,
                const double *Y, const int incY,
                const double *beta, const double *W, const int incW,
                const double *Z, const int incZ,
                double *A, const int lda);
void ATL_zger2c(const int M, const int N,
                const double *alpha, const double *X, const int incX,
                const double *Y, const int incY,
                const double *beta, const double *W, const int incW,
                const double *Z, const int incZ,
                double *A, const int lda);
void ATL_zher(const enum ATLAS_UPLO Uplo, const int N, const double alpha,
              const double *X, const int incX, double *A, const int lda);
void ATL_zhpr(const enum ATLAS_UPLO Uplo, const int N, const double alpha,
                   const double *X, const int incX, double *A);
void ATL_zher2(const enum ATLAS_UPLO Uplo, const int N,
               const double *alpha, const double *X, const int incX,
               const double *Y, const int incY, double *A, const int lda);
void ATL_zhpr2(const enum ATLAS_UPLO Uplo, const int N,
               const double *alpha, const double *X, const int incX,
               const double *Y, const int incY, double *Ap);


#endif
