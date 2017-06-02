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
#ifndef ATLAS_REFLEVEL3_H
#define ATLAS_REFLEVEL3_H

#include "atlas_enum.h"
/*
 * =====================================================================
 * Prototypes for Level 3 Reference ATLAS BLAS routines
 * =====================================================================
 */
void       ATL_srefgemm
(
  const enum ATLAS_TRANS, const enum ATLAS_TRANS,
  const int,              const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsymm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyrk
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const float,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyr2k
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_sreftrmm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_drefgemm
(
  const enum ATLAS_TRANS, const enum ATLAS_TRANS,
  const int,              const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsymm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyrk
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const double,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyr2k
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_dreftrmm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_crefgemm
(
  const enum ATLAS_TRANS, const enum ATLAS_TRANS,
  const int,              const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhemm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefherk
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const float,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_crefher2k
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_crefsymm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsyrk
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsyr2k
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_creftrmm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_zrefgemm
(
  const enum ATLAS_TRANS, const enum ATLAS_TRANS,
  const int,              const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhemm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefherk
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const double,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_zrefher2k
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_zrefsymm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsyrk
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsyr2k
(
  const enum ATLAS_UPLO,  const enum ATLAS_TRANS,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zreftrmm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsm
(
  const enum ATLAS_SIDE,  const enum ATLAS_UPLO,
  const enum ATLAS_TRANS, const enum ATLAS_DIAG,
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

#endif
/*
 * End of atlas_reflevel3.h
 */
