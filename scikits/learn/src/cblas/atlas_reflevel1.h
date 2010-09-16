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

#ifndef ATLAS_REFLEVEL1_H
#define ATLAS_REFLEVEL1_H

/*
 * =====================================================================
 * Prototypes for Level 1 Reference ATLAS BLAS routines
 * =====================================================================
 */
void       ATL_srefrotg
(
   float *,
   float *,
   float *,
   float *
);

void       ATL_srefrotmg
(
   float *,
   float *,
   float *,
   const float,
   float *
);

float      ATL_srefnrm2
(
   const int,
   const float *,          const int
);

float      ATL_srefasum
(
   const int,
   const float *,          const int
);

int        ATL_isrefamax
(
   const int,
   const float *,          const int
);

void       ATL_srefscal
(
   const int,
   const float,
   float *,                const int
);

void       ATL_srefswap
(
   const int,
   float *,                const int,
   float *,                const int
);

void       ATL_srefcopy
(
   const int,
   const float *,          const int,
   float *,                const int
);

void       ATL_srefaxpy
(
   const int,
   const float,
   const float *,          const int,
   float *,                const int
);

void       ATL_srefrot
(
   const int,
   float *,                const int,
   float *,                const int,
   const float,
   const float
);

void       ATL_srefrotm
(
   const int,
   float *,                const int,
   float *,                const int,
   const float *
);

float      ATL_srefdot
(
   const int,
   const float *,          const int,
   const float *,          const int
);

float      ATL_sdsrefdot
(
   const int,
   const float,
   const float *,          const int,
   const float *,          const int
);

double     ATL_dsrefdot
(
   const int,
   const float *,          const int,
   const float *,          const int
);

void       ATL_drefrotg
(
   double *,
   double *,
   double *,
   double *
);

void       ATL_drefrotmg
(
   double *,
   double *,
   double *,
   const double,
   double *
);

double     ATL_drefnrm2
(
   const int,
   const double *,         const int
);

double     ATL_drefasum
(
   const int,
   const double *,         const int
);

int        ATL_idrefamax
(
   const int,
   const double *,         const int
);

void       ATL_drefscal
(
   const int,
   const double,
   double *,               const int
);

void       ATL_drefswap
(
   const int,
   double *,               const int,
   double *,               const int
);

void       ATL_drefcopy
(
   const int,
   const double *,         const int,
   double *,               const int
);

void       ATL_drefaxpy
(
   const int,
   const double,
   const double *,         const int,
   double *,               const int
);

void       ATL_drefrot
(
   const int,
   double *,               const int,
   double *,               const int,
   const double,
   const double
);

void       ATL_drefrotm
(
   const int,
   double *,               const int,
   double *,               const int,
   const double *
);

double     ATL_drefdot
(
   const int,
   const double *,         const int,
   const double *,         const int
);

void       ATL_crefrotg
(
   float *,
   const float *,
   float *,
   float *
);

float      ATL_screfnrm2
(
   const int,
   const float *,          const int
);

float      ATL_screfasum
(
   const int,
   const float *,          const int
);

int        ATL_icrefamax
(
   const int,
   const float *,          const int
);

void       ATL_crefscal
(
   const int,
   const float *,
   float *,                const int
);

void       ATL_csrefscal
(
   const int,
   const float,
   float *,                const int
);

void       ATL_crefswap
(
   const int,
   float *,                const int,
   float *,                const int
);

void       ATL_crefcopy
(
   const int,
   const float *,          const int,
   float *,                const int
);

void       ATL_crefaxpy
(
   const int,
   const float *,
   const float *,          const int,
   float *,                const int
);

void       ATL_csrefrot
(
   const int,
   float *,                const int,
   float *,                const int,
   const float,
   const float
);

void       ATL_crefdotc_sub
(
   const int,
   const float *,          const int,
   const float *,          const int,
   float *
);

void       ATL_crefdotu_sub
(
   const int,
   const float *,          const int,
   const float *,          const int,
   float *
);

void       ATL_zrefrotg
(
   double *,
   const double *,
   double *,
   double *
);

double     ATL_dzrefnrm2
(
   const int,
   const double *,         const int
);

double     ATL_dzrefasum
(
   const int,
   const double *,         const int
);

int        ATL_izrefamax
(
   const int,
   const double *,         const int
);

void       ATL_zrefscal
(
   const int,
   const double *,
   double *,               const int
);

void       ATL_zdrefscal
(
   const int,
   const double,
   double *,               const int
);

void       ATL_zrefswap
(
   const int,
   double *,               const int,
   double *,               const int
);

void       ATL_zrefcopy
(
   const int,
   const double *,         const int,
   double *,               const int
);

void       ATL_zrefaxpy
(
   const int,
   const double *,
   const double *,         const int,
   double *,               const int
);

void       ATL_zdrefrot
(
   const int,
   double *,               const int,
   double *,               const int,
   const double,
   const double
);

void       ATL_zrefdotc_sub
(
   const int,
   const double *,         const int,
   const double *,         const int,
   double *
);

void       ATL_zrefdotu_sub
(
   const int,
   const double *,         const int,
   const double *,         const int,
   double *
);

#endif
/*
 * End of atlas_reflevel1.h
 */
