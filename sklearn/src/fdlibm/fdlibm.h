
/* @(#)fdlibm.h 1.5 04/04/22 */
/*
 * ====================================================
 * Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice 
 * is preserved.
 * ====================================================
 */

/* Sometimes it's necessary to define __LITTLE_ENDIAN explicitly
   but these catch some common cases. */

#if defined(i386) || defined(i486) || \
	defined(intel) || defined(x86) || defined(i86pc) || \
	defined(__alpha) || defined(__osf__)
#define __LITTLE_ENDIAN
#endif

#ifdef __LITTLE_ENDIAN
#define __HI(x) *(1+(int*)&x)
#define __LO(x) *(int*)&x
#define __HIp(x) *(1+(int*)x)
#define __LOp(x) *(int*)x
#else
#define __HI(x) *(int*)&x
#define __LO(x) *(1+(int*)&x)
#define __HIp(x) *(int*)x
#define __LOp(x) *(1+(int*)x)
#endif

#ifndef _CDEFS_H_
#ifdef __STDC__
#define	__P(p)	p
#else
#define	__P(p)	()
#endif
#endif

/*
 * ANSI/POSIX
 */

extern int signgam;

#ifndef MAXFLOAT
#define	MAXFLOAT	((float)3.40282346638528860e+38)
#endif


enum fdversion {fdlibm_ieee = -1, fdlibm_svid, fdlibm_xopen, fdlibm_posix};

#define _LIB_VERSION_TYPE enum fdversion
#define _LIB_VERSION _fdlib_version  

/* if global variable _LIB_VERSION is not desirable, one may 
 * change the following to be a constant by: 
 *	#define _LIB_VERSION_TYPE const enum version
 * In that case, after one initializes the value _LIB_VERSION (see
 * s_lib_version.c) during compile time, it cannot be modified
 * in the middle of a program
 */ 
extern  _LIB_VERSION_TYPE  _LIB_VERSION;

#define _IEEE_  fdlibm_ieee
#define _SVID_  fdlibm_svid
#define _XOPEN_ fdlibm_xopen
#define _POSIX_ fdlibm_posix

struct fdlibm_exception {
	int type;
	char *name;
	double arg1;
	double arg2;
	double retval;
};

#define	HUGE		MAXFLOAT

/* 
 * set X_TLOSS = pi*2**52, which is possibly defined in <values.h>
 * (one may replace the following line by "#include <values.h>")
 */

#define X_TLOSS		1.41484755040568800000e+16 

#ifndef DOMAIN
#define	DOMAIN		1
#endif
#ifndef SING
#define	SING		2
#endif
#ifndef OVERFLOW
#define	OVERFLOW	3
#endif
#ifndef UNDERFLOW
#define	UNDERFLOW	4
#endif
#ifndef TLOSS
#define	TLOSS		5
#endif
#ifndef PLOSS
#define	PLOSS		6
#endif

/*
 * ANSI/POSIX
 */
/*
extern double acos __P((double));
extern double asin __P((double));
extern double atan __P((double));
extern double atan2 __P((double, double));
extern double cos __P((double));
extern double sin __P((double));
extern double tan __P((double));

extern double cosh __P((double));
extern double sinh __P((double));
extern double tanh __P((double));

extern double exp __P((double));
extern double frexp __P((double, int *));
extern double ldexp __P((double, int));
extern double log __P((double));
extern double log10 __P((double));
extern double modf __P((double, double *));

extern double pow __P((double, double));
extern double sqrt __P((double));

extern double ceil __P((double));
*/
extern double fabs __P((double));
extern double floor __P((double));
/*
extern double fmod __P((double, double));

extern double erf __P((double));
extern double erfc __P((double));
extern double gamma __P((double));
extern double hypot __P((double, double));
extern int isnan __P((double));
*/
extern int finite __P((double));
/*
extern double j0 __P((double));
extern double j1 __P((double));
extern double jn __P((int, double));
*/
extern double lgamma __P((double));
/*
extern double y0 __P((double));
extern double y1 __P((double));
extern double yn __P((int, double));

extern double acosh __P((double));
extern double asinh __P((double));
extern double atanh __P((double));
extern double cbrt __P((double));
extern double logb __P((double));
extern double nextafter __P((double, double));
extern double remainder __P((double, double));
#ifdef _SCALB_INT
extern double scalb __P((double, int));
#else
extern double scalb __P((double, double));
#endif

extern int matherr __P((struct exception *));

*/
/*
 * IEEE Test Vector
 */

/*
extern double significand __P((double));
*/

/*
 * Functions callable from C, intended to support IEEE arithmetic.
 */


extern double copysign __P((double, double));
/*
extern int ilogb __P((double));
*/
extern double rint __P((double));
/*
extern double scalbn __P((double, int));
*/

/*
 * BSD math library entry points
 */
/*
extern double expm1 __P((double));
extern double log1p __P((double));
*/
/*
 * Reentrant version of gamma & lgamma; passes signgam back by reference
 * as the second argument; user must allocate space for signgam.
 */
#ifdef _REENTRANT
extern double gamma_r __P((double, int *));
extern double lgamma_r __P((double, int *));
#endif	/* _REENTRANT */

/* ieee style elementary functions */
/*extern double __ieee754_sqrt __P((double));			
extern double __ieee754_acos __P((double));			
extern double __ieee754_acosh __P((double));			
*/
extern double __ieee754_log __P((double));			
/*extern double __ieee754_atanh __P((double));			
extern double __ieee754_asin __P((double));			
extern double __ieee754_atan2 __P((double,double));			
extern double __ieee754_exp __P((double));
extern double __ieee754_cosh __P((double));
extern double __ieee754_fmod __P((double,double));
extern double __ieee754_pow __P((double,double));
*/
extern double __ieee754_lgamma_r __P((double,int *));
/*extern double __ieee754_gamma_r __P((double,int *));*/
extern double __ieee754_lgamma __P((double));
/*
extern double __ieee754_gamma __P((double));
extern double __ieee754_log10 __P((double));
extern double __ieee754_sinh __P((double));
extern double __ieee754_hypot __P((double,double));
extern double __ieee754_j0 __P((double));
extern double __ieee754_j1 __P((double));
extern double __ieee754_y0 __P((double));
extern double __ieee754_y1 __P((double));
extern double __ieee754_jn __P((int,double));
extern double __ieee754_yn __P((int,double));
extern double __ieee754_remainder __P((double,double));
extern int    __ieee754_rem_pio2 __P((double,double*));
#ifdef _SCALB_INT
extern double __ieee754_scalb __P((double,int));
#else
extern double __ieee754_scalb __P((double,double));
#endif
*/
/* fdlibm kernel function */
extern double __kernel_standard __P((double,double,int));	
extern double __kernel_sin __P((double,double,int));
extern double __kernel_cos __P((double,double));
/*
extern double __kernel_tan __P((double,double,int));
extern int    __kernel_rem_pio2 __P((double*,double*,int,int,int,const int*));*/

