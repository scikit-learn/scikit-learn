/*
 *             Automatically Tuned Linear Algebra Software v3.6.0
 *                    (C) Copyright 1997 R. Clint Whaley
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
#ifndef ATLAS_ENUM_H
   #define ATLAS_ENUM_H

   #define CBLAS_ENUM_ONLY
   #include "cblas.h"
   #undef CBLAS_ENUM_ONLY

   #define ATLAS_ORDER CBLAS_ORDER
      #define AtlasRowMajor CblasRowMajor
      #define AtlasColMajor CblasColMajor
   #define ATLAS_TRANS CBLAS_TRANSPOSE
      #define AtlasNoTrans CblasNoTrans
      #define AtlasTrans CblasTrans
      #define AtlasConjTrans CblasConjTrans
   #define ATLAS_UPLO CBLAS_UPLO
      #define AtlasUpper CblasUpper
      #define AtlasLower CblasLower
   #define ATLAS_DIAG CBLAS_DIAG
      #define AtlasNonUnit CblasNonUnit
      #define AtlasUnit CblasUnit
   #define ATLAS_SIDE CBLAS_SIDE
      #define AtlasLeft  CblasLeft
      #define AtlasRight CblasRight

#endif

