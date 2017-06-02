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
#ifndef ATLAS_REFLVL3_H
#define ATLAS_REFLVL3_H
/*
 * =====================================================================
 * Prototypes for Level 3 Reference Internal ATLAS BLAS routines
 * =====================================================================
 */
void       ATL_srefgemmNN
(
  const int,              const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgemmNT
(
  const int,              const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgemmTN
(
  const int,              const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefgemmTT
(
  const int,              const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsymmLL
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsymmLU
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsymmRL
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsymmRU
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyrkLN
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyrkLT
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyrkUN
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyrkUT
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyr2kLN
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyr2kLT
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyr2kUN
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_srefsyr2kUT
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_sreftrmmLLNN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmLLNU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmLLTN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmLLTU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmLUNN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmLUNU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmLUTN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmLUTU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmRLNN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmRLNU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmRLTN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmRLTU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmRUNN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmRUNU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmRUTN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrmmRUTU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmLLNN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmLLNU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmLLTN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmLLTU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmLUNN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmLUNU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmLUTN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmLUTU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmRLNN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmRLNU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmRLTN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmRLTU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmRUNN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmRUNU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmRUTN
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_sreftrsmRUTU
(
  const int,              const int,
  const float,
  const float *,          const int,
  float *,                const int
);

void       ATL_drefgemmNN
(
  const int,              const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgemmNT
(
  const int,              const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgemmTN
(
  const int,              const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefgemmTT
(
  const int,              const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsymmLL
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsymmLU
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsymmRL
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsymmRU
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyrkLN
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyrkLT
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyrkUN
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyrkUT
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyr2kLN
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyr2kLT
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyr2kUN
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_drefsyr2kUT
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_dreftrmmLLNN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmLLNU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmLLTN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmLLTU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmLUNN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmLUNU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmLUTN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmLUTU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmRLNN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmRLNU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmRLTN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmRLTU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmRUNN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmRUNU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmRUTN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrmmRUTU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmLLNN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmLLNU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmLLTN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmLLTU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmLUNN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmLUNU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmLUTN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmLUTU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmRLNN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmRLNU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmRLTN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmRLTU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmRUNN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmRUNU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmRUTN
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_dreftrsmRUTU
(
  const int,              const int,
  const double,
  const double *,         const int,
  double *,               const int
);

void       ATL_crefgemmNN
(
  const int,              const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemmNT
(
  const int,              const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemmNC
(
  const int,              const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemmTN
(
  const int,              const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemmTT
(
  const int,              const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemmTC
(
  const int,              const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemmCN
(
  const int,              const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemmCT
(
  const int,              const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefgemmCC
(
  const int,              const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhemmLL
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhemmLU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhemmRL
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefhemmRU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefherkLN
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_crefherkLC
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_crefherkUN
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_crefherkUC
(
  const int,              const int,
  const float,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_crefher2kLN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_crefher2kLC
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_crefher2kUN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_crefher2kUC
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float,
  float *,                const int
);

void       ATL_crefsymmLL
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsymmLU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsymmRL
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsymmRU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsyrkLN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsyrkLT
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsyrkUN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsyrkUT
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsyr2kLN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsyr2kLT
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsyr2kUN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_crefsyr2kUT
(
  const int,              const int,
  const float *,
  const float *,          const int,
  const float *,          const int,
  const float *,
  float *,                const int
);

void       ATL_creftrmmLLNN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmLLNU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmLLTN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmLLTU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmLLCN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmLLCU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmLUNN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmLUNU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmLUTN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmLUTU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmLUCN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmLUCU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmRLNN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmRLNU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmRLTN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmRLTU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmRLCN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmRLCU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmRUNN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmRUNU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmRUTN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmRUTU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmRUCN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrmmRUCU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmLLNN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmLLNU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmLLTN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmLLTU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmLLCN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmLLCU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmLUNN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmLUNU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmLUTN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmLUTU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmLUCN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmLUCU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmRLNN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmRLNU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmRLTN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmRLTU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmRLCN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmRLCU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmRUNN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmRUNU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmRUTN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmRUTU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmRUCN
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_creftrsmRUCU
(
  const int,              const int,
  const float *,
  const float *,          const int,
  float *,                const int
);

void       ATL_zrefgemmNN
(
  const int,              const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemmNT
(
  const int,              const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemmNC
(
  const int,              const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemmTN
(
  const int,              const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemmTT
(
  const int,              const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemmTC
(
  const int,              const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemmCN
(
  const int,              const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemmCT
(
  const int,              const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefgemmCC
(
  const int,              const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhemmLL
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhemmLU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhemmRL
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefhemmRU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefherkLN
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_zrefherkLC
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_zrefherkUN
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_zrefherkUC
(
  const int,              const int,
  const double,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_zrefher2kLN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_zrefher2kLC
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_zrefher2kUN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_zrefher2kUC
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double,
  double *,               const int
);

void       ATL_zrefsymmLL
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsymmLU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsymmRL
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsymmRU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsyrkLN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsyrkLT
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsyrkUN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsyrkUT
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsyr2kLN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsyr2kLT
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsyr2kUN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zrefsyr2kUT
(
  const int,              const int,
  const double *,
  const double *,         const int,
  const double *,         const int,
  const double *,
  double *,               const int
);

void       ATL_zreftrmmLLNN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmLLNU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmLLTN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmLLTU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmLLCN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmLLCU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmLUNN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmLUNU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmLUTN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmLUTU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmLUCN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmLUCU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmRLNN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmRLNU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmRLTN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmRLTU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmRLCN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmRLCU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmRUNN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmRUNU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmRUTN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmRUTU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmRUCN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrmmRUCU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmLLNN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmLLNU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmLLTN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmLLTU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmLLCN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmLLCU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmLUNN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmLUNU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmLUTN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmLUTU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmLUCN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmLUCU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmRLNN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmRLNU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmRLTN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmRLTU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmRLCN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmRLCU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmRUNN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmRUNU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmRUTN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmRUTU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmRUCN
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

void       ATL_zreftrsmRUCU
(
  const int,              const int,
  const double *,
  const double *,         const int,
  double *,               const int
);

#endif
/*
 * End of atlas_reflvl3.h
 */
