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
#ifndef ATLAS_REFLVL2_H
#define ATLAS_REFLVL2_H
/*
 * =====================================================================
 * Prototypes for Level 2 Reference Internal ATLAS BLAS routines
 * =====================================================================
 */
void       ATL_srefgbmvN
(
  const int,              const int,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgbmvT
(
  const int,              const int,
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgpmvUN
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgpmvUT
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgpmvLN
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgpmvLT
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgemvN
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgemvT
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgprL
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefgprU
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefsbmvL
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsbmvU
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefspmvL
(
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefspmvU
(
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsprL
(
  const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefsprU
(
  const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefspr2L
(
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefspr2U
(
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefsymvL
(
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsymvU
(
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyrL
(
  const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefsyrU
(
  const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefsyr2L
(
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_srefsyr2U
(
  const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbmvLNN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbmvLNU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbmvLTN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbmvLTU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbmvUNN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbmvUNU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbmvUTN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbmvUTU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpmvLNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpmvLNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpmvLTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpmvLTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpmvUNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpmvUNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpmvUTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpmvUTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmvLNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmvLNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmvLTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmvLTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmvUNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmvUNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmvUTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmvUTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbsvLNN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbsvLNU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbsvLTN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbsvLTU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbsvUNN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbsvUNU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbsvUTN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftbsvUTU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpsvLNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpsvLNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpsvLTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpsvLTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpsvUNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpsvUNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpsvUTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftpsvUTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsvLNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsvLNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsvLTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsvLTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsvUNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsvUNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsvUTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsvUTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_drefgbmvN
(
  const int,              const int,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgbmvT
(
  const int,              const int,
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgpmvUN
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgpmvUT
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgpmvLN
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgpmvLT
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgemvN
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgemvT
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgprL
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefgprU
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefsbmvL
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsbmvU
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefspmvL
(
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefspmvU
(
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsprL
(
  const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefsprU
(
  const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefspr2L
(
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefspr2U
(
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefsymvL
(
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsymvU
(
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyrL
(
  const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefsyrU
(
  const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefsyr2L
(
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_drefsyr2U
(
  const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbmvLNN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbmvLNU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbmvLTN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbmvLTU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbmvUNN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbmvUNU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbmvUTN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbmvUTU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpmvLNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpmvLNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpmvLTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpmvLTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpmvUNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpmvUNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpmvUTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpmvUTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmvLNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmvLNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmvLTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmvLTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmvUNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmvUNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmvUTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmvUTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbsvLNN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbsvLNU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbsvLTN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbsvLTU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbsvUNN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbsvUNU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbsvUTN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftbsvUTU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpsvLNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpsvLNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpsvLTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpsvLTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpsvUNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpsvUNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpsvUTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftpsvUTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsvLNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsvLNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsvLTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsvLTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsvUNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsvUNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsvUTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsvUTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_crefgbmvN
(
  const int,              const int,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgbmvT
(
  const int,              const int,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgbmvC
(
  const int,              const int,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgbmvH
(
  const int,              const int,
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgpmvUN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgpmvUT
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgpmvUC
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgpmvUH
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgpmvLN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgpmvLT
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgpmvLC
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgpmvLH
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemvN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemvT
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemvC
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemvH
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgprcL
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefgprcU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefgpruL
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefgpruU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefhbmvL
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhbmvU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhpmvL
(
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhpmvU
(
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhprL
(
  const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefhprU
(
  const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefhpr2L
(
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefhpr2U
(
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefhemvL
(
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhemvU
(
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefherL
(
  const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefherU
(
  const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefher2L
(
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_crefher2U
(
  const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvLNN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvLNU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvLTN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvLTU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvLCN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvLCU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvLHN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvLHU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvUNN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvUNU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvUTN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvUTU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvUCN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvUCU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvUHN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbmvUHU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvLNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvLNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvLTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvLTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvLCN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvLCU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvLHN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvLHU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvUNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvUNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvUTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvUTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvUCN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvUCU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvUHN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpmvUHU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvLNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvLNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvLTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvLTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvLCN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvLCU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvLHN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvLHU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvUNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvUNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvUTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvUTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvUCN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvUCU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvUHN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmvUHU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvLNN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvLNU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvLTN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvLTU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvLCN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvLCU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvLHN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvLHU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvUNN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvUNU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvUTN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvUTU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvUCN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvUCU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvUHN
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftbsvUHU
(
  const int,              const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvLNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvLNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvLTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvLTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvLCN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvLCU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvLHN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvLHU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvUNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvUNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvUTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvUTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvUCN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvUCU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvUHN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftpsvUHU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvLNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvLNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvLTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvLTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvLCN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvLCU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvLHN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvLHU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvUNN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvUNU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvUTN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvUTU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvUCN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvUCU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvUHN
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsvUHU
(
  const int,
  const float *,          const int,
  float *,                const int
);

void       ATL_zrefgbmvN
(
  const int,              const int,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgbmvT
(
  const int,              const int,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgbmvC
(
  const int,              const int,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgbmvH
(
  const int,              const int,
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgpmvUN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgpmvUT
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgpmvUC
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgpmvUH
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgpmvLN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgpmvLT
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgpmvLC
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgpmvLH
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemvN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemvT
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemvC
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemvH
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgprcL
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefgprcU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefgpruL
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefgpruU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefhbmvL
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhbmvU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhpmvL
(
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhpmvU
(
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhprL
(
  const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefhprU
(
  const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefhpr2L
(
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefhpr2U
(
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefhemvL
(
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhemvU
(
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefherL
(
  const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefherU
(
  const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefher2L
(
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zrefher2U
(
  const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvLNN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvLNU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvLTN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvLTU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvLCN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvLCU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvLHN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvLHU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvUNN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvUNU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvUTN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvUTU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvUCN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvUCU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvUHN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbmvUHU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvLNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvLNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvLTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvLTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvLCN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvLCU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvLHN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvLHU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvUNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvUNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvUTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvUTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvUCN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvUCU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvUHN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpmvUHU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvLNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvLNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvLTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvLTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvLCN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvLCU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvLHN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvLHU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvUNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvUNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvUTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvUTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvUCN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvUCU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvUHN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmvUHU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvLNN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvLNU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvLTN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvLTU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvLCN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvLCU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvLHN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvLHU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvUNN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvUNU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvUTN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvUTU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvUCN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvUCU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvUHN
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftbsvUHU
(
  const int,              const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvLNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvLNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvLTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvLTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvLCN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvLCU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvLHN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvLHU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvUNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvUNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvUTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvUTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvUCN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvUCU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvUHN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftpsvUHU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvLNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvLNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvLTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvLTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvLCN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvLCU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvLHN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvLHU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvUNN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvUNU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvUTN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvUTU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvUCN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvUCU
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvUHN
(
  const int,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsvUHU
(
  const int,
  const double *,         const int,
  double *,               const int
);

#endif
/*
 * End of atlas_reflvl2.h
 */
