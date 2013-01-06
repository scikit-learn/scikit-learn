
/* @(#)k_standard.c 1.3 95/01/18 */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice 
 * is preserved.
 * ====================================================
 *
 */

#include "fdlibm.h"
#include <errno.h>

#ifndef _USE_WRITE
#include <stdio.h>			/* fputs(), stderr */
#define	WRITE2(u,v)	fputs(u, stderr)
#else	/* !defined(_USE_WRITE) */
#include <unistd.h>			/* write */
#define	WRITE2(u,v)	write(2, u, v)
#undef fflush
#endif	/* !defined(_USE_WRITE) */

static double zero = 0.0;	/* used as const */

/* 
 * Standard conformance (non-IEEE) on exception cases.
 * Mapping:
 *	1 -- acos(|x|>1)
 *	2 -- asin(|x|>1)
 *	3 -- atan2(+-0,+-0)
 *	4 -- hypot overflow
 *	5 -- cosh overflow
 *	6 -- exp overflow
 *	7 -- exp underflow
 *	8 -- y0(0)
 *	9 -- y0(-ve)
 *	10-- y1(0)
 *	11-- y1(-ve)
 *	12-- yn(0)
 *	13-- yn(-ve)
 *	14-- lgamma(finite) overflow
 *	15-- lgamma(-integer)
 *	16-- log(0)
 *	17-- log(x<0)
 *	18-- log10(0)
 *	19-- log10(x<0)
 *	20-- pow(0.0,0.0)
 *	21-- pow(x,y) overflow
 *	22-- pow(x,y) underflow
 *	23-- pow(0,negative) 
 *	24-- pow(neg,non-integral)
 *	25-- sinh(finite) overflow
 *	26-- sqrt(negative)
 *      27-- fmod(x,0)
 *      28-- remainder(x,0)
 *	29-- acosh(x<1)
 *	30-- atanh(|x|>1)
 *	31-- atanh(|x|=1)
 *	32-- scalb overflow
 *	33-- scalb underflow
 *	34-- j0(|x|>X_TLOSS)
 *	35-- y0(x>X_TLOSS)
 *	36-- j1(|x|>X_TLOSS)
 *	37-- y1(x>X_TLOSS)
 *	38-- jn(|x|>X_TLOSS, n)
 *	39-- yn(x>X_TLOSS, n)
 *	40-- gamma(finite) overflow
 *	41-- gamma(-integer)
 *	42-- pow(NaN,0.0)
 */


#ifdef __STDC__
	double __kernel_standard(double x, double y, int type) 
#else
	double __kernel_standard(x,y,type) 
	double x,y; int type;
#endif
{
	struct fdlibm_exception exc;
#ifndef HUGE_VAL	/* this is the only routine that uses HUGE_VAL */ 
#define HUGE_VAL inf
	double inf = 0.0;

	__HI(inf) = 0x7ff00000;	/* set inf to infinite */
#endif

#ifdef _USE_WRITE
	(void) fflush(stdout);
#endif
	exc.arg1 = x;
	exc.arg2 = y;
	switch(type) {
	    case 1:
		/* acos(|x|>1) */
		exc.type = DOMAIN;
		exc.name = "acos";
		exc.retval = zero;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if(_LIB_VERSION == _SVID_) {
		    (void) WRITE2("acos: DOMAIN error\n", 19);
		  }
		  errno = EDOM;
		}
		break;
	    case 2:
		/* asin(|x|>1) */
		exc.type = DOMAIN;
		exc.name = "asin";
		exc.retval = zero;
		if(_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if(_LIB_VERSION == _SVID_) {
		    	(void) WRITE2("asin: DOMAIN error\n", 19);
		  }
		  errno = EDOM;
		}
		break;
	    case 3:
		/* atan2(+-0,+-0) */
		exc.arg1 = y;
		exc.arg2 = x;
		exc.type = DOMAIN;
		exc.name = "atan2";
		exc.retval = zero;
		if(_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if(_LIB_VERSION == _SVID_) {
			(void) WRITE2("atan2: DOMAIN error\n", 20);
		      }
		  errno = EDOM;
		}
		break;
	    case 4:
		/* hypot(finite,finite) overflow */
		exc.type = OVERFLOW;
		exc.name = "hypot";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = HUGE;
		else
		  exc.retval = HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = ERANGE;
		else if (!matherr(&exc)) {
			errno = ERANGE;
		}
		break;
	    case 5:
		/* cosh(finite) overflow */
		exc.type = OVERFLOW;
		exc.name = "cosh";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = HUGE;
		else
		  exc.retval = HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = ERANGE;
		else if (!matherr(&exc)) {
			errno = ERANGE;
		}
		break;
	    case 6:
		/* exp(finite) overflow */
		exc.type = OVERFLOW;
		exc.name = "exp";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = HUGE;
		else
		  exc.retval = HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = ERANGE;
		else if (!matherr(&exc)) {
			errno = ERANGE;
		}
		break;
	    case 7:
		/* exp(finite) underflow */
		exc.type = UNDERFLOW;
		exc.name = "exp";
		exc.retval = zero;
		if (_LIB_VERSION == _POSIX_)
		  errno = ERANGE;
		else if (!matherr(&exc)) {
			errno = ERANGE;
		}
		break;
	    case 8:
		/* y0(0) = -inf */
		exc.type = DOMAIN;	/* should be SING for IEEE */
		exc.name = "y0";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = -HUGE;
		else
		  exc.retval = -HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("y0: DOMAIN error\n", 17);
		      }
		  errno = EDOM;
		}
		break;
	    case 9:
		/* y0(x<0) = NaN */
		exc.type = DOMAIN;
		exc.name = "y0";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = -HUGE;
		else
		  exc.retval = -HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("y0: DOMAIN error\n", 17);
		      }
		  errno = EDOM;
		}
		break;
	    case 10:
		/* y1(0) = -inf */
		exc.type = DOMAIN;	/* should be SING for IEEE */
		exc.name = "y1";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = -HUGE;
		else
		  exc.retval = -HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("y1: DOMAIN error\n", 17);
		      }
		  errno = EDOM;
		}
		break;
	    case 11:
		/* y1(x<0) = NaN */
		exc.type = DOMAIN;
		exc.name = "y1";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = -HUGE;
		else
		  exc.retval = -HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("y1: DOMAIN error\n", 17);
		      }
		  errno = EDOM;
		}
		break;
	    case 12:
		/* yn(n,0) = -inf */
		exc.type = DOMAIN;	/* should be SING for IEEE */
		exc.name = "yn";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = -HUGE;
		else
		  exc.retval = -HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("yn: DOMAIN error\n", 17);
		      }
		  errno = EDOM;
		}
		break;
	    case 13:
		/* yn(x<0) = NaN */
		exc.type = DOMAIN;
		exc.name = "yn";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = -HUGE;
		else
		  exc.retval = -HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("yn: DOMAIN error\n", 17);
		      }
		  errno = EDOM;
		}
		break;
	    case 14:
		/* lgamma(finite) overflow */
		exc.type = OVERFLOW;
		exc.name = "lgamma";
                if (_LIB_VERSION == _SVID_)
                  exc.retval = HUGE;
                else
                  exc.retval = HUGE_VAL;
                if (_LIB_VERSION == _POSIX_)
			errno = ERANGE;
                else if (!matherr(&exc)) {
                        errno = ERANGE;
		}
		break;
	    case 15:
		/* lgamma(-integer) or lgamma(0) */
		exc.type = SING;
		exc.name = "lgamma";
                if (_LIB_VERSION == _SVID_)
                  exc.retval = HUGE;
                else
                  exc.retval = HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("lgamma: SING error\n", 19);
		      }
		  errno = EDOM;
		}
		break;
	    case 16:
		/* log(0) */
		exc.type = SING;
		exc.name = "log";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = -HUGE;
		else
		  exc.retval = -HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = ERANGE;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("log: SING error\n", 16);
		      }
		  errno = EDOM;
		}
		break;
	    case 17:
		/* log(x<0) */
		exc.type = DOMAIN;
		exc.name = "log";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = -HUGE;
		else
		  exc.retval = -HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("log: DOMAIN error\n", 18);
		      }
		  errno = EDOM;
		}
		break;
	    case 18:
		/* log10(0) */
		exc.type = SING;
		exc.name = "log10";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = -HUGE;
		else
		  exc.retval = -HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = ERANGE;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("log10: SING error\n", 18);
		      }
		  errno = EDOM;
		}
		break;
	    case 19:
		/* log10(x<0) */
		exc.type = DOMAIN;
		exc.name = "log10";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = -HUGE;
		else
		  exc.retval = -HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("log10: DOMAIN error\n", 20);
		      }
		  errno = EDOM;
		}
		break;
	    case 20:
		/* pow(0.0,0.0) */
		/* error only if _LIB_VERSION == _SVID_ */
		exc.type = DOMAIN;
		exc.name = "pow";
		exc.retval = zero;
		if (_LIB_VERSION != _SVID_) exc.retval = 1.0;
		else if (!matherr(&exc)) {
			(void) WRITE2("pow(0,0): DOMAIN error\n", 23);
			errno = EDOM;
		}
		break;
	    case 21:
		/* pow(x,y) overflow */
		exc.type = OVERFLOW;
		exc.name = "pow";
		if (_LIB_VERSION == _SVID_) {
		  exc.retval = HUGE;
		  y *= 0.5;
		  if(x<zero&&rint(y)!=y) exc.retval = -HUGE;
		} else {
		  exc.retval = HUGE_VAL;
		  y *= 0.5;
		  if(x<zero&&rint(y)!=y) exc.retval = -HUGE_VAL;
		}
		if (_LIB_VERSION == _POSIX_)
		  errno = ERANGE;
		else if (!matherr(&exc)) {
			errno = ERANGE;
		}
		break;
	    case 22:
		/* pow(x,y) underflow */
		exc.type = UNDERFLOW;
		exc.name = "pow";
		exc.retval =  zero;
		if (_LIB_VERSION == _POSIX_)
		  errno = ERANGE;
		else if (!matherr(&exc)) {
			errno = ERANGE;
		}
		break;
	    case 23:
		/* 0**neg */
		exc.type = DOMAIN;
		exc.name = "pow";
		if (_LIB_VERSION == _SVID_) 
		  exc.retval = zero;
		else
		  exc.retval = -HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("pow(0,neg): DOMAIN error\n", 25);
		      }
		  errno = EDOM;
		}
		break;
	    case 24:
		/* neg**non-integral */
		exc.type = DOMAIN;
		exc.name = "pow";
		if (_LIB_VERSION == _SVID_) 
		    exc.retval = zero;
		else 
		    exc.retval = zero/zero;	/* X/Open allow NaN */
		if (_LIB_VERSION == _POSIX_) 
		   errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("neg**non-integral: DOMAIN error\n", 32);
		      }
		  errno = EDOM;
		}
		break;
	    case 25:
		/* sinh(finite) overflow */
		exc.type = OVERFLOW;
		exc.name = "sinh";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = ( (x>zero) ? HUGE : -HUGE);
		else
		  exc.retval = ( (x>zero) ? HUGE_VAL : -HUGE_VAL);
		if (_LIB_VERSION == _POSIX_)
		  errno = ERANGE;
		else if (!matherr(&exc)) {
			errno = ERANGE;
		}
		break;
	    case 26:
		/* sqrt(x<0) */
		exc.type = DOMAIN;
		exc.name = "sqrt";
		if (_LIB_VERSION == _SVID_)
		  exc.retval = zero;
		else
		  exc.retval = zero/zero;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("sqrt: DOMAIN error\n", 19);
		      }
		  errno = EDOM;
		}
		break;
            case 27:
                /* fmod(x,0) */
                exc.type = DOMAIN;
                exc.name = "fmod";
                if (_LIB_VERSION == _SVID_)
                    exc.retval = x;
		else
		    exc.retval = zero/zero;
                if (_LIB_VERSION == _POSIX_)
                  errno = EDOM;
                else if (!matherr(&exc)) {
                  if (_LIB_VERSION == _SVID_) {
                    (void) WRITE2("fmod:  DOMAIN error\n", 20);
                  }
                  errno = EDOM;
                }
                break;
            case 28:
                /* remainder(x,0) */
                exc.type = DOMAIN;
                exc.name = "remainder";
                exc.retval = zero/zero;
                if (_LIB_VERSION == _POSIX_)
                  errno = EDOM;
                else if (!matherr(&exc)) {
                  if (_LIB_VERSION == _SVID_) {
                    (void) WRITE2("remainder: DOMAIN error\n", 24);
                  }
                  errno = EDOM;
                }
                break;
            case 29:
                /* acosh(x<1) */
                exc.type = DOMAIN;
                exc.name = "acosh";
                exc.retval = zero/zero;
                if (_LIB_VERSION == _POSIX_)
                  errno = EDOM;
                else if (!matherr(&exc)) {
                  if (_LIB_VERSION == _SVID_) {
                    (void) WRITE2("acosh: DOMAIN error\n", 20);
                  }
                  errno = EDOM;
                }
                break;
            case 30:
                /* atanh(|x|>1) */
                exc.type = DOMAIN;
                exc.name = "atanh";
                exc.retval = zero/zero;
                if (_LIB_VERSION == _POSIX_)
                  errno = EDOM;
                else if (!matherr(&exc)) {
                  if (_LIB_VERSION == _SVID_) {
                    (void) WRITE2("atanh: DOMAIN error\n", 20);
                  }
                  errno = EDOM;
                }
                break;
            case 31:
                /* atanh(|x|=1) */
                exc.type = SING;
                exc.name = "atanh";
		exc.retval = x/zero;	/* sign(x)*inf */
                if (_LIB_VERSION == _POSIX_)
                  errno = EDOM;
                else if (!matherr(&exc)) {
                  if (_LIB_VERSION == _SVID_) {
                    (void) WRITE2("atanh: SING error\n", 18);
                  }
                  errno = EDOM;
                }
                break;
	    case 32:
		/* scalb overflow; SVID also returns +-HUGE_VAL */
		exc.type = OVERFLOW;
		exc.name = "scalb";
		exc.retval = x > zero ? HUGE_VAL : -HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = ERANGE;
		else if (!matherr(&exc)) {
			errno = ERANGE;
		}
		break;
	    case 33:
		/* scalb underflow */
		exc.type = UNDERFLOW;
		exc.name = "scalb";
		exc.retval = copysign(zero,x);
		if (_LIB_VERSION == _POSIX_)
		  errno = ERANGE;
		else if (!matherr(&exc)) {
			errno = ERANGE;
		}
		break;
	    case 34:
		/* j0(|x|>X_TLOSS) */
                exc.type = TLOSS;
                exc.name = "j0";
                exc.retval = zero;
                if (_LIB_VERSION == _POSIX_)
                        errno = ERANGE;
                else if (!matherr(&exc)) {
                        if (_LIB_VERSION == _SVID_) {
                                (void) WRITE2(exc.name, 2);
                                (void) WRITE2(": TLOSS error\n", 14);
                        }
                        errno = ERANGE;
                }        
		break;
	    case 35:
		/* y0(x>X_TLOSS) */
                exc.type = TLOSS;
                exc.name = "y0";
                exc.retval = zero;
                if (_LIB_VERSION == _POSIX_)
                        errno = ERANGE;
                else if (!matherr(&exc)) {
                        if (_LIB_VERSION == _SVID_) {
                                (void) WRITE2(exc.name, 2);
                                (void) WRITE2(": TLOSS error\n", 14);
                        }
                        errno = ERANGE;
                }        
		break;
	    case 36:
		/* j1(|x|>X_TLOSS) */
                exc.type = TLOSS;
                exc.name = "j1";
                exc.retval = zero;
                if (_LIB_VERSION == _POSIX_)
                        errno = ERANGE;
                else if (!matherr(&exc)) {
                        if (_LIB_VERSION == _SVID_) {
                                (void) WRITE2(exc.name, 2);
                                (void) WRITE2(": TLOSS error\n", 14);
                        }
                        errno = ERANGE;
                }        
		break;
	    case 37:
		/* y1(x>X_TLOSS) */
                exc.type = TLOSS;
                exc.name = "y1";
                exc.retval = zero;
                if (_LIB_VERSION == _POSIX_)
                        errno = ERANGE;
                else if (!matherr(&exc)) {
                        if (_LIB_VERSION == _SVID_) {
                                (void) WRITE2(exc.name, 2);
                                (void) WRITE2(": TLOSS error\n", 14);
                        }
                        errno = ERANGE;
                }        
		break;
	    case 38:
		/* jn(|x|>X_TLOSS) */
                exc.type = TLOSS;
                exc.name = "jn";
                exc.retval = zero;
                if (_LIB_VERSION == _POSIX_)
                        errno = ERANGE;
                else if (!matherr(&exc)) {
                        if (_LIB_VERSION == _SVID_) {
                                (void) WRITE2(exc.name, 2);
                                (void) WRITE2(": TLOSS error\n", 14);
                        }
                        errno = ERANGE;
                }        
		break;
	    case 39:
		/* yn(x>X_TLOSS) */
                exc.type = TLOSS;
                exc.name = "yn";
                exc.retval = zero;
                if (_LIB_VERSION == _POSIX_)
                        errno = ERANGE;
                else if (!matherr(&exc)) {
                        if (_LIB_VERSION == _SVID_) {
                                (void) WRITE2(exc.name, 2);
                                (void) WRITE2(": TLOSS error\n", 14);
                        }
                        errno = ERANGE;
                }        
		break;
	    case 40:
		/* gamma(finite) overflow */
		exc.type = OVERFLOW;
		exc.name = "gamma";
                if (_LIB_VERSION == _SVID_)
                  exc.retval = HUGE;
                else
                  exc.retval = HUGE_VAL;
                if (_LIB_VERSION == _POSIX_)
		  errno = ERANGE;
                else if (!matherr(&exc)) {
                  errno = ERANGE;
                }
		break;
	    case 41:
		/* gamma(-integer) or gamma(0) */
		exc.type = SING;
		exc.name = "gamma";
                if (_LIB_VERSION == _SVID_)
                  exc.retval = HUGE;
                else
                  exc.retval = HUGE_VAL;
		if (_LIB_VERSION == _POSIX_)
		  errno = EDOM;
		else if (!matherr(&exc)) {
		  if (_LIB_VERSION == _SVID_) {
			(void) WRITE2("gamma: SING error\n", 18);
		      }
		  errno = EDOM;
		}
		break;
	    case 42:
		/* pow(NaN,0.0) */
		/* error only if _LIB_VERSION == _SVID_ & _XOPEN_ */
		exc.type = DOMAIN;
		exc.name = "pow";
		exc.retval = x;
		if (_LIB_VERSION == _IEEE_ ||
		    _LIB_VERSION == _POSIX_) exc.retval = 1.0;
		else if (!matherr(&exc)) {
			errno = EDOM;
		}
		break;
	}
	return exc.retval; 
}
