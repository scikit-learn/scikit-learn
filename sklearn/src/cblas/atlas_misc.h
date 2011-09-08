/*
 *             Automatically Tuned Linear Algebra Software v3.9.25
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

#include <stdio.h>
#include <stdlib.h>
#include "atlas_enum.h"

#ifndef ATLAS_MISC_H
#define ATLAS_MISC_H
#include "atlas_type.h"
#ifdef ATL_PROFILE
   extern int ATL_ProfGemmCameFrom;
#endif
/*
 * If using a C99 compiler, define the restrict attribute, which says that
 * this is the only pointer referencing the given section of memory
 */
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif

/*
 * Some useful macro functions
 */
#if (defined(PentiumCPS) || defined(ATL_USEPTHREADS)) && !defined(WALL)
   #define WALL
#endif
#ifndef time00
   #if defined(WALL)
      #define time00 ATL_walltime
   #else
      #define time00 ATL_cputime
   #endif
#endif
#define Mabs(x) ( (x) >= 0 ? (x) : -(x) )
#define Mmax(x, y) ( (x) > (y) ? (x) : (y) )
#define Mmin(x, y) ( (x) > (y) ? (y) : (x) )
#define Mlowcase(C) ( ((C) > 64 && (C) < 91) ? (C) | 32 : (C) )
#define Mupcase(C) ( ((C) > 96 && (C) < 123) ? (C) & 0xDF : (C) )
/*
 * packed indexing functions (upper & lower)
 */

#define Mjoin(pre, nam) my_join(pre, nam)
#define my_join(pre, nam) pre ## nam
#define Mstr2(m) # m
#define Mstr(m) Mstr2(m)

#define ATL_assert(n_) \
{ \
   if (!(n_)) \
   { \
      ATL_xerbla(0, __FILE__, "assertion %s failed, line %d of file %s\n", \
                 Mstr(n_), __LINE__, __FILE__); \
   } \
}

/*
 * Define some C99 features that we use when we know the compiler supports them
 */
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define INLINE inline
   #define RESTRICT restrict
#else
   #define INLINE
   #define RESTRICT
#endif

#ifdef ATL_LONG_INT
   #define ATL_INT long
   #define ATL_CINT const long
#else
   #define ATL_INT int
   #define ATL_CINT const int
#endif


#define ATL_QTYPE long double
#if defined(SREAL)
   #define EPS 5.0e-7
   #define TYPE float
   #define PRE s
   #define UPR s
   #define PREU S
   #define PATL ATL_s
   #define PATU ATLU_s
   #define UATL ATLU_s
   #define CBLA cblas_s
   #define PATLU ATL_s
   #define ATL_rone   1.0f
   #define ATL_rnone -1.0f
   #define ATL_rzero   0.0f
   #define ATL_typify(m_) Mjoin(m_,f)
   #include "atlas_ssysinfo.h"
#elif defined(DREAL)
   #define EPS 1.0e-15
   #define TYPE double
   #define PRE d
   #define UPR d
   #define PREU D
   #define PATL ATL_d
   #define PATU ATLU_d
   #define UATL ATLU_d
   #define CBLA cblas_d
   #define PATLU ATL_d
   #define ATL_rone   1.0
   #define ATL_rnone -1.0
   #define ATL_rzero   0.0
   #define ATL_typify(m_) m_
   #include "atlas_dsysinfo.h"
#elif defined (QREAL)
   #define EPS 1.9259299443872358530559779425849273E-34L
   #define TYPE long double
   #define PRE q
   #define UPR q
   #define PREU Q
   #define PATL ATL_q
   #define PATU ATLU_q
   #define CBLA cblas_q
   #define ATL_rone 1.0
   #define ATL_rnone -1.0
   #define ATL_rzero 0.0
#elif defined(SCPLX)
   #define EPS 5.0e-7
   #define TYPE float
   #define PRE c
   #define UPR s
   #define PREU C
   #define PATL ATL_c
   #define PATLU ATL_s
   #define PATU ATLU_c
   #define UATL ATLU_s
   #define ATL_rone  1.0f
   #define ATL_rnone -1.0f
   #define ATL_rzero   0.0f
   #define ATL_typify(m_) Mjoin(m_,f)
   #define CBLA cblas_c
   #include "atlas_csysinfo.h"
#elif defined(DCPLX)
   #define TYPE double
   #define PRE z
   #define UPR d
   #define PREU Z
   #define PATL ATL_z
   #define PATLU ATL_d
   #define PATU ATLU_z
   #define UATL ATLU_d
   #define EPS 1.0e-15
   #define ATL_rone   1.0
   #define ATL_rnone -1.0
   #define ATL_rzero   0.0
   #define ATL_typify(m_) m_
   #define CBLA cblas_z
   #include "atlas_zsysinfo.h"
#elif defined(QCPLX)
   #define TYPE ATL_QTYPE
   #define PRE e
   #deffine UPR q
   #define PREU E
   #define PATL ATL_e
   #define PATLU ATL_q
   #define UATL ATLU_q
   #define EPS 1.9259299443872358530559779425849273E-34L
   #define ATL_rone 1.0
   #define ATL_rnone -1.0
   #define ATL_rzero 0.0
   #define CBLA cblas_e
#endif

#if defined (SREAL) || defined (DREAL) || defined (SCPLX) || defined (DCPLX)
   #define ATL_sizeof Mjoin(PATL,size)
   #define ATL_MulBySize Mjoin(PATL,MulBySize)
   #define ATL_DivBySize Mjoin(PATL,DivBySize)
#else
   #define ATL_sizeof sizeof(TYPE)
   #define ATL_MulBySize Mjoin(Mjoin(ATL_,PRE),MulBySize)(N_) \
      ((N_)*sizeof(TYPE))
   #define ATL_DivBySize Mjoin(Mjoin(ATL_,PRE),DivBySize)(N_) \
      ((N_)/sizeof(TYPE))
#endif

#if ( defined(SREAL) || defined(DREAL) || defined(QREAL) )
   #define TREAL
   #define SHIFT
   #define SCALAR TYPE
   #define SADD &
   #define SVAL
   #define SVVAL *
   #define SCALAR_IS_ONE(M_scalar) ((M_scalar) == ATL_rone)
   #define SCALAR_IS_NONE(M_scalar) ((M_scalar) == ATL_rnone)
   #define SCALAR_IS_ZERO(M_scalar) ((M_scalar) == ATL_rzero)
#elif defined(SCPLX) || defined(DCPLX) || defined(QCPLX)
   #define TCPLX
/*
 * c = b*c + v;
 */
   #define CMULT2(v, a, b, tmp) \
   { \
      tmp = *(a) * *(b) - *(a+1) * *(b+1); \
      *(b+1) = *(a) * *(b+1) + *(a+1) * *(b) + *(v+1); \
      *(b) = tmp + *v; \
   }
   #define SHIFT << 1
   #define SCALAR TYPE *
   #define SADD
   #define SVAL *
   #define SVVAL
   #define SCALAR_IS_ONE(M_scalar) \
      ( (*(M_scalar) == ATL_rone) && ((M_scalar)[1] == ATL_rzero) )
   #define SCALAR_IS_NONE(M_scalar) \
      ( (*(M_scalar) == ATL_rnone) && ((M_scalar)[1] == ATL_rzero) )
   #define SCALAR_IS_ZERO(M_scalar) \
      ( (*(M_scalar) == ATL_rzero) && ((M_scalar)[1] == ATL_rzero) )
#endif

#if defined(ALPHA1)
   #define ATL_MulByALPHA(x_) (x_)
   #define NM _a1
#elif defined (ALPHA0)
   #define ATL_MulByALPHA(x_) ATL_rzero
   #define NM _a0
#elif defined (ALPHAN1)
   #define ATL_MulByALPHA(x_) (-(x_))
   #define NM _an1
#elif defined (ALPHAXI0)
   #define ATL_MulByALPHA(x_) (ralpha*(x_))
   #define NM _aXi0
#elif defined (ALPHA1C)
   #define NM _a1c
#elif defined (ALPHAN1C)
   #define NM _an1c
#elif defined (ALPHAXI0C)
   #define NM _aXi0c
#elif defined (ALPHAXC)
   #define NM _aXc
#elif defined (ALPHAX)
   #define ATL_MulByALPHA(x_) (alpha*(x_))
   #define NM _aX
#endif

#if defined(BETA1)
   #define ATL_MulByBETA(x_) (x_)
   #define MSTAT A[i] += v[i]
   #define BNM _b1
#elif defined(BETA1C)
   #define BNM _b1c
#elif defined(BETAN1)
   #define ATL_MulByBETA(x_) (-(x_))
   #define MSTAT A[i] = v[i] - A[i]
   #define BNM _bn1
#elif defined(BETAN1C)
   #define BNM _bn1c
#elif defined(BETA0)
   #define ATL_MulByBETA(x_) ATL_rzero
   #define MSTAT A[i] = v[i]
   #define BNM _b0
#elif defined (BETAXI0)
   #define BNM _bXi0
   #define ATL_MulByBETA(x_) (rbeta*(x_))
#elif defined (BETAXI0C)
   #define BNM _bXi0c
#elif defined (BETAX)
   #define ATL_MulByBETA(x_) (beta*(x_))
   #define MSTAT A[i] = beta*A[i] + v[i]
   #define BNM _bX
#elif defined (BETAXC)
   #define BNM _bXc
#endif

/* any alignment below this forces data copy in gemm */
#ifndef ATL_MinMMAlign
   #define ATL_MinMMAlign 16
#endif
#if (ATL_MinMMAlign == 1 || ATL_MinMMAlign == 0)
   #define ATL_DataIsMinAligned(ptr) 1
#elif (ATL_MinMMAlign == 2)
   #define ATL_DataIsMinAligned(ptr) \
      ( (((size_t) (ptr))>>1)<<1 == (size_t) (ptr) )
#elif (ATL_MinMMAlign == 4)
   #define ATL_DataIsMinAligned(ptr) \
      ( (((size_t) (ptr))>>2)<<2 == (size_t) (ptr) )
#elif (ATL_MinMMAlign == 8)
   #define ATL_DataIsMinAligned(ptr) \
      ( (((size_t) (ptr))>>3)<<3 == (size_t) (ptr) )
#elif (ATL_MinMMAlign == 16)
   #define ATL_DataIsMinAligned(ptr) \
      ( (((size_t) (ptr))>>4)<<4 == (size_t) (ptr) )
#elif (ATL_MinMMAlign == 32)
   #define ATL_DataIsMinAligned(ptr) \
      ( (((size_t) (ptr))>>5)<<5 == (size_t) (ptr) )
#elif (ATL_MinMMAlign == 64)
   #define ATL_DataIsMinAligned(ptr) \
      ( (((size_t) (ptr))>>6)<<6 == (size_t) (ptr) )
#elif (ATL_MinMMAlign == 128)
   #define ATL_DataIsMinAligned(ptr) \
      ( (((size_t) (ptr))>>7)<<7 == (size_t) (ptr) )
#else
   #define ATL_DataIsMinAligned(ptr) \
      ( (((size_t) (ptr))/ATL_MinMMAlign)*ATL_MinMMAlign == (size_t) (ptr) )
#endif

#define ATL_Cachelen 32
#if (ATL_Cachelen == 4)
   #define ATL_MulByCachelen(N_) ( (N_) << 2 )
   #define ATL_DivByCachelen(N_) ( (N_) >> 2 )
#elif (ATL_Cachelen == 8)
   #define ATL_MulByCachelen(N_) ( (N_) << 3 )
   #define ATL_DivByCachelen(N_) ( (N_) >> 3 )
#elif (ATL_Cachelen == 16)
   #define ATL_MulByCachelen(N_) ( (N_) << 4 )
   #define ATL_DivByCachelen(N_) ( (N_) >> 4 )
#elif (ATL_Cachelen == 32)
   #define ATL_MulByCachelen(N_) ( (N_) << 5 )
   #define ATL_DivByCachelen(N_) ( (N_) >> 5 )
#elif (ATL_Cachelen == 64)
   #define ATL_MulByCachelen(N_) ( (N_) << 6 )
   #define ATL_DivByCachelen(N_) ( (N_) >> 6 )
#elif (ATL_Cachelen == 128)
   #define ATL_MulByCachelen(N_) ( (N_) << 7 )
   #define ATL_DivByCachelen(N_) ( (N_) >> 7 )
#elif (ATL_Cachelen == 256)
   #define ATL_MulByCachelen(N_) ( (N_) << 8 )
   #define ATL_DivByCachelen(N_) ( (N_) >> 8 )
#else
   #define ATL_MulByCachelen(N_) ( (N_) * ATL_Cachelen )
   #define ATL_DivByCachelen(N_) ( (N_) / ATL_Cachelen )
#endif

#if (ATL_Cachelen < ATL_MinMMAlign)
   Force a compilation error if our required alignment is at least the
   minimum!!@^
#endif

#define ATL_AlignPtr(vp) \
   (void*) (ATL_Cachelen + ATL_MulByCachelen(ATL_DivByCachelen((size_t) (vp))))

#define ATL_FindPtrAdjust(vp, iadj_) \
{ \
   (iadj_) = ((size_t)(vp))-ATL_MulByCachelen(ATL_DivByCachelen((size_t)(vp)));\
   if (iadj_) \
   { \
      if ( (iadj_) == ATL_MulBySize(ATL_DivBySize(iadj_)) ) \
         (iadj_) = ATL_DivBySize(iadj_); \
      else (iadj_) = 0; \
   }\
}
#define ATL_FindMatAdjust(vp_, lda_, iadj_) \
{ \
   if (ATL_MulByCachelen(ATL_DivByCachelen(ATL_MulBySize(lda_))) \
       == ATL_MulBySize(lda_)) \
   { \
      ATL_FindPtrAdjust(vp_, iadj_); \
   } \
   else (iadj_) = 0; \
}

#define ATL_sqrtLL(x, res) \
   asm ("fsqrt" : "=t" (res) : "0" (x));

/*
 * Find N necessary for alignment.  Written as function for optimization,
 * declared static to encourage inlining
 */
static int ATL_AlignOffset
(const int N,       /* max return value */
 const void *vp,    /* pointer to be aligned */
 const int inc,     /* size of each elt, in bytes */
 const int align)   /* required alignment, in bytes */
{
   const int p = align/inc;
   const size_t k=(size_t)vp, j=k/inc;
   int iret;
   if (k == (j)*inc && p*inc == align)
   {
      iret = ((j+p-1) / p)*p - j;
      if (iret <= N) return(iret);
   }
   return(N);
}
static void *ATL_Align2Ptr(const void *pu, const void *pA)
/*
 * Aligns pu%ATL_Cachelen to pA%ATL_Cachelen by adding at most ATL_Cachlen
 * to pu
 * RETURNS: pu possibly incremented so that it has same alignment as pA
 */
{
   size_t tu = (size_t) pu, ta = (size_t) pA;

   tu -= (tu/ATL_Cachelen)*ATL_Cachelen;
   ta -= (ta/ATL_Cachelen)*ATL_Cachelen;
   if (tu <= ta)
      tu = ta;
   else
      tu = ta + ATL_Cachelen-tu;
   tu += (size_t) pu;
   return((void*) tu);
}


/*
 * Gcc links in crap that MSVC++ and DVF can't handle if you use stdout
 * or stderr, so use this beautiful kludge to avoid this problem -- RCW
 */
#ifdef GCCWIN

#include <stdarg.h>
static int WINFPRINTF(FILE *fpout, char *form, ...)
{
   int ierr=0;
   va_list argptr;

   va_start(argptr, form);
   if (fpout == NULL) ierr = vprintf(form, argptr);
   else ierr = vfprintf(fpout, form, argptr);
   va_end(argptr);

   return(ierr);
}

#ifdef stdout
   #undef stdout
#endif
#ifdef stderr
   #undef stderr
#endif
#ifdef assert
   #undef assert
#endif

#define stdout NULL
#define stderr NULL
#define fprintf WINFPRINTF
#define assert WINASSERT
#define WINASSERT(n_) \
{ \
   if (!(n_)) \
   { \
      printf("assertion %s failed, line %d of file %s\n", \
             Mstr(n_), __LINE__, __FILE__); \
      exit(1); \
   } \
}

#endif

#include "atlas_aux.h"

#endif
