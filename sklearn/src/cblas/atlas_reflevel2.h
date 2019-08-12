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

#ifndef ATLAS_REFLEVEL2_H
#define ATLAS_REFLEVEL2_H

#include "atlas_enum.h"
/*
 * =====================================================================
 * Prototypes for Level 2 Reference ATLAS BLAS routines
 * =====================================================================
 */
void       ATL_srefgbmv
(
  const ATLAS_TRANS,
  const int,              const int,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgpmv
(
  const ATLAS_UPLO,
  const ATLAS_TRANS,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgemv
(
  const ATLAS_TRANS,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgpr
(
  const ATLAS_UPLO,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefger
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefsbmv
(
  const ATLAS_UPLO,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefspmv
(
  const ATLAS_UPLO,
  const int,
  const float,
  const float *,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefspr
(
  const ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  float *
);

void       ATL_srefspr2
(
  const ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *
);

void       ATL_srefsymv
(
  const ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyr
(
  const ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefsyr2
(
  const ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbmv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbsv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpmv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const float *,
  float *,                const int
);

void       ATL_sreftpsv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const float *,
  float *,                const int
);

void       ATL_sreftrmv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_drefgbmv
(
  const ATLAS_TRANS,
  const int,              const int,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgpmv
(
  const ATLAS_UPLO,
  const ATLAS_TRANS,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgemv
(
  const ATLAS_TRANS,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgpr
(
  const ATLAS_UPLO,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefger
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefsbmv
(
  const ATLAS_UPLO,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefspmv
(
  const ATLAS_UPLO,
  const int,
  const double,
  const double *,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefspr
(
  const ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  double *
);

void       ATL_drefspr2
(
  const ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *
);

void       ATL_drefsymv
(
  const ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyr
(
  const ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefsyr2
(
  const ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbmv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbsv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpmv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const double *,
  double *,               const int
);

void       ATL_dreftpsv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const double *,
  double *,               const int
);

void       ATL_dreftrmv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_crefgbmv
(
  const ATLAS_TRANS,
  const int,              const int,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgpmv
(
  const ATLAS_UPLO,
  const ATLAS_TRANS,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemv
(
  const ATLAS_TRANS,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgprc
(
  const ATLAS_UPLO,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefgpru
(
  const ATLAS_UPLO,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefgerc
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefgeru
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefhbmv
(
  const ATLAS_UPLO,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhpmv
(
  const ATLAS_UPLO,
  const int,
  const float *,
  const float *,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhpr
(
  const ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  float *
);

void       ATL_crefhpr2
(
  const ATLAS_UPLO,
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *
);

void       ATL_crefhemv
(
  const ATLAS_UPLO,
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefher
(
  const ATLAS_UPLO,
  const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefher2
(
  const ATLAS_UPLO,
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const float *,
  float *,                const int
);

void       ATL_creftpsv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const float *,
  float *,                const int
);

void       ATL_creftrmv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_zrefgbmv
(
  const ATLAS_TRANS,
  const int,              const int,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgpmv
(
  const ATLAS_UPLO,
  const ATLAS_TRANS,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemv
(
  const ATLAS_TRANS,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgprc
(
  const ATLAS_UPLO,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefgpru
(
  const ATLAS_UPLO,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefgerc
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefgeru
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefhbmv
(
  const ATLAS_UPLO,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhpmv
(
  const ATLAS_UPLO,
  const int,
  const double *,
  const double *,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhpr
(
  const ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  double *
);

void       ATL_zrefhpr2
(
  const ATLAS_UPLO,
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *
);

void       ATL_zrefhemv
(
  const ATLAS_UPLO,
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefher
(
  const ATLAS_UPLO,
  const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefher2
(
  const ATLAS_UPLO,
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const double *,
  double *,               const int
);

void       ATL_zreftpsv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const double *,
  double *,               const int
);

void       ATL_zreftrmv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsv
(
  const ATLAS_UPLO,  const ATLAS_TRANS, const ATLAS_DIAG,
  const int,
  const double *,         const int,
  double *,               const int
);

#endif
/*
 * End of atlas_reflevel2.h
 */
