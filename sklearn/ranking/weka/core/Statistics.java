package weka.core;

/**
 * Class implementing some distributions, tests, etc. The code is mostly adapted from the CERN
 * Jet Java libraries:
 * 
 * Copyright 2001 University of Waikato
 * Copyright 1999 CERN - European Organization for Nuclear Research.
 * Permission to use, copy, modify, distribute and sell this software and its documentation for
 * any purpose is hereby granted without fee, provided that the above copyright notice appear
 * in all copies and that both that copyright notice and this permission notice appear in
 * supporting documentation. 
 * CERN and the University of Waikato make no representations about the suitability of this 
 * software for any purpose. It is provided "as is" without expressed or implied warranty.
 *
 * @author peter.gedeck@pharma.Novartis.com
 * @author wolfgang.hoschek@cern.ch
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @version $Revision: 5616 $
 */
public class Statistics
  implements RevisionHandler {

  /** Some constants */
  protected static final double MACHEP =  1.11022302462515654042E-16;
  protected static final double MAXLOG =  7.09782712893383996732E2;
  protected static final double MINLOG = -7.451332191019412076235E2;
  protected static final double MAXGAM = 171.624376956302725;
  protected static final double SQTPI  =  2.50662827463100050242E0;
  protected static final double SQRTH  =  7.07106781186547524401E-1;
  protected static final double LOGPI  =  1.14472988584940017414;
  
  protected static final double big    =  4.503599627370496e15;
  protected static final double biginv =  2.22044604925031308085e-16;

  /*************************************************
   *    COEFFICIENTS FOR METHOD  normalInverse()   *
   *************************************************/
  /* approximation for 0 <= |y - 0.5| <= 3/8 */
  protected static final double P0[] = {
    -5.99633501014107895267E1,
    9.80010754185999661536E1,
    -5.66762857469070293439E1,
    1.39312609387279679503E1,
    -1.23916583867381258016E0,
  };
  protected static final double Q0[] = {
    /* 1.00000000000000000000E0,*/
    1.95448858338141759834E0,
    4.67627912898881538453E0,
    8.63602421390890590575E1,
    -2.25462687854119370527E2,
    2.00260212380060660359E2,
    -8.20372256168333339912E1,
    1.59056225126211695515E1,
    -1.18331621121330003142E0,
  };
  
  /* Approximation for interval z = sqrt(-2 log y ) between 2 and 8
   * i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
   */
  protected static final double P1[] = {
    4.05544892305962419923E0,
    3.15251094599893866154E1,
    5.71628192246421288162E1,
    4.40805073893200834700E1,
    1.46849561928858024014E1,
    2.18663306850790267539E0,
    -1.40256079171354495875E-1,
    -3.50424626827848203418E-2,
    -8.57456785154685413611E-4,
  };
  protected static final double Q1[] = {
    /*  1.00000000000000000000E0,*/
    1.57799883256466749731E1,
    4.53907635128879210584E1,
    4.13172038254672030440E1,
    1.50425385692907503408E1,
    2.50464946208309415979E0,
    -1.42182922854787788574E-1,
    -3.80806407691578277194E-2,
    -9.33259480895457427372E-4,
  };
  
  /* Approximation for interval z = sqrt(-2 log y ) between 8 and 64
   * i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
   */
  protected static final double  P2[] = {
    3.23774891776946035970E0,
    6.91522889068984211695E0,
    3.93881025292474443415E0,
    1.33303460815807542389E0,
    2.01485389549179081538E-1,
    1.23716634817820021358E-2,
    3.01581553508235416007E-4,
    2.65806974686737550832E-6,
    6.23974539184983293730E-9,
  };
  protected static final double  Q2[] = {
    /*  1.00000000000000000000E0,*/
    6.02427039364742014255E0,
    3.67983563856160859403E0,
    1.37702099489081330271E0,
    2.16236993594496635890E-1,
    1.34204006088543189037E-2,
    3.28014464682127739104E-4,
    2.89247864745380683936E-6,
    6.79019408009981274425E-9,
  };
  
  /**
   * Computes standard error for observed values of a binomial
   * random variable.
   *
   * @param p the probability of success
   * @param n the size of the sample
   * @return the standard error
   */
  public static double binomialStandardError(double p, int n) {
    
    if (n == 0) {
      return 0; 
    }
    return Math.sqrt((p*(1-p))/(double) n);
  }
  
  /**
   * Returns chi-squared probability for given value and degrees
   * of freedom. (The probability that the chi-squared variate
   * will be greater than x for the given degrees of freedom.)
   *
   * @param x the value
   * @param v the number of degrees of freedom
   * @return the chi-squared probability
   */
  public static double chiSquaredProbability(double x, double v) { 

    if( x < 0.0 || v < 1.0 ) return 0.0;
    return incompleteGammaComplement( v/2.0, x/2.0 );
  }

  /**
   * Computes probability of F-ratio.
   *
   * @param F the F-ratio
   * @param df1 the first number of degrees of freedom
   * @param df2 the second number of degrees of freedom
   * @return the probability of the F-ratio.
   */
  public static double FProbability(double F, int df1, int df2) {
    
    return incompleteBeta( df2/2.0, df1/2.0, df2/(df2+df1*F) );
  }

  /**
   * Returns the area under the Normal (Gaussian) probability density
   * function, integrated from minus infinity to <tt>x</tt>
   * (assumes mean is zero, variance is one).
   * <pre>
   *                            x
   *                             -
   *                   1        | |          2
   *  normal(x)  = ---------    |    exp( - t /2 ) dt
   *               sqrt(2pi)  | |
   *                           -
   *                          -inf.
   *
   *             =  ( 1 + erf(z) ) / 2
   *             =  erfc(z) / 2
   * </pre>
   * where <tt>z = x/sqrt(2)</tt>.
   * Computation is via the functions <tt>errorFunction</tt> and <tt>errorFunctionComplement</tt>.
   *
   * @param a the z-value
   * @return the probability of the z value according to the normal pdf
   */
  public static double normalProbability(double a) { 

    double x, y, z;
 
    x = a * SQRTH;
    z = Math.abs(x);
 
    if( z < SQRTH ) y = 0.5 + 0.5 * errorFunction(x);
    else {
      y = 0.5 * errorFunctionComplemented(z);
      if( x > 0 )  y = 1.0 - y;
    } 
    return y;
  }

  /**
   * Returns the value, <tt>x</tt>, for which the area under the
   * Normal (Gaussian) probability density function (integrated from
   * minus infinity to <tt>x</tt>) is equal to the argument <tt>y</tt>
   * (assumes mean is zero, variance is one).
   * <p>
   * For small arguments <tt>0 < y < exp(-2)</tt>, the program computes
   * <tt>z = sqrt( -2.0 * log(y) )</tt>;  then the approximation is
   * <tt>x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z)</tt>.
   * There are two rational functions P/Q, one for <tt>0 < y < exp(-32)</tt>
   * and the other for <tt>y</tt> up to <tt>exp(-2)</tt>. 
   * For larger arguments,
   * <tt>w = y - 0.5</tt>, and  <tt>x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2))</tt>.
   *
   * @param y0 the area under the normal pdf
   * @return the z-value
   */
  public static double normalInverse(double y0) { 

    double x, y, z, y2, x0, x1;
    int code;

    final double s2pi = Math.sqrt(2.0*Math.PI);

    if( y0 <= 0.0 ) throw new IllegalArgumentException();
    if( y0 >= 1.0 ) throw new IllegalArgumentException();
    code = 1;
    y = y0;
    if( y > (1.0 - 0.13533528323661269189) ) { /* 0.135... = exp(-2) */
      y = 1.0 - y;
      code = 0;
    }

    if( y > 0.13533528323661269189 ) {
      y = y - 0.5;
      y2 = y * y;
      x = y + y * (y2 * polevl( y2, P0, 4)/p1evl( y2, Q0, 8 ));
      x = x * s2pi; 
      return(x);
    }

    x = Math.sqrt( -2.0 * Math.log(y) );
    x0 = x - Math.log(x)/x;

    z = 1.0/x;
    if( x < 8.0 ) /* y > exp(-32) = 1.2664165549e-14 */
      x1 = z * polevl( z, P1, 8 )/p1evl( z, Q1, 8 );
    else
      x1 = z * polevl( z, P2, 8 )/p1evl( z, Q2, 8 );
    x = x0 - x1;
    if( code != 0 )
      x = -x;
    return( x );
  }

  /**
   * Returns natural logarithm of gamma function.
   *
   * @param x the value
   * @return natural logarithm of gamma function
   */
  public static double lnGamma(double x) {

    double p, q, w, z;
 
    double A[] = {
      8.11614167470508450300E-4,
      -5.95061904284301438324E-4,
      7.93650340457716943945E-4,
      -2.77777777730099687205E-3,
      8.33333333333331927722E-2
    };
    double B[] = {
      -1.37825152569120859100E3,
      -3.88016315134637840924E4,
      -3.31612992738871184744E5,
      -1.16237097492762307383E6,
      -1.72173700820839662146E6,
      -8.53555664245765465627E5
    };
    double C[] = {
      /* 1.00000000000000000000E0, */
      -3.51815701436523470549E2,
      -1.70642106651881159223E4,
      -2.20528590553854454839E5,
      -1.13933444367982507207E6,
      -2.53252307177582951285E6,
      -2.01889141433532773231E6
    };
 
    if( x < -34.0 ) {
      q = -x;
      w = lnGamma(q);
      p = Math.floor(q);
      if( p == q ) throw new ArithmeticException("lnGamma: Overflow");
      z = q - p;
      if( z > 0.5 ) {
	p += 1.0;
	z = p - q;
      }
      z = q * Math.sin( Math.PI * z );
      if( z == 0.0 ) throw new 
	ArithmeticException("lnGamma: Overflow");
      z = LOGPI - Math.log( z ) - w;
      return z;
    }
 
    if( x < 13.0 ) {
      z = 1.0;
      while( x >= 3.0 ) {
	x -= 1.0;
	z *= x;
      }
      while( x < 2.0 ) {
	if( x == 0.0 ) throw new 
	  ArithmeticException("lnGamma: Overflow");
	z /= x;
	x += 1.0;
      }
      if( z < 0.0 ) z = -z;
      if( x == 2.0 ) return Math.log(z);
      x -= 2.0;
      p = x * polevl( x, B, 5 ) / p1evl( x, C, 6);
      return( Math.log(z) + p );
    }
 
    if( x > 2.556348e305 ) throw new ArithmeticException("lnGamma: Overflow");
 
    q = ( x - 0.5 ) * Math.log(x) - x + 0.91893853320467274178;
  
    if( x > 1.0e8 ) return( q );
 
    p = 1.0/(x*x);
    if( x >= 1000.0 )
      q += ((   7.9365079365079365079365e-4 * p
		- 2.7777777777777777777778e-3) *p
	    + 0.0833333333333333333333) / x;
    else
      q += polevl( p, A, 4 ) / x;
    return q;
  }

  /**
   * Returns the error function of the normal distribution.
   * The integral is
   * <pre>
   *                           x 
   *                            -
   *                 2         | |          2
   *   erf(x)  =  --------     |    exp( - t  ) dt.
   *              sqrt(pi)   | |
   *                          -
   *                           0
   * </pre>
   * <b>Implementation:</b>
   * For <tt>0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2)</tt>; otherwise
   * <tt>erf(x) = 1 - erfc(x)</tt>.
   * <p>
   * Code adapted from the <A HREF="http://www.sci.usq.edu.au/staff/leighb/graph/Top.html">
   * Java 2D Graph Package 2.4</A>,
   * which in turn is a port from the
   * <A HREF="http://people.ne.mediaone.net/moshier/index.html#Cephes">Cephes 2.2</A>
   * Math Library (C).
   *
   * @param a the argument to the function.
   */
  public static double errorFunction(double x) { 
    double y, z;
    final double T[] = {
      9.60497373987051638749E0,
      9.00260197203842689217E1,
      2.23200534594684319226E3,
      7.00332514112805075473E3,
      5.55923013010394962768E4
    };
    final double U[] = {
      //1.00000000000000000000E0,
      3.35617141647503099647E1,
      5.21357949780152679795E2,
      4.59432382970980127987E3,
      2.26290000613890934246E4,
      4.92673942608635921086E4
    };
  
    if( Math.abs(x) > 1.0 ) return( 1.0 - errorFunctionComplemented(x) );
    z = x * x;
    y = x * polevl( z, T, 4 ) / p1evl( z, U, 5 );
    return y;
  }

  /**
   * Returns the complementary Error function of the normal distribution.
   * <pre>
   *  1 - erf(x) =
   *
   *                           inf. 
   *                             -
   *                  2         | |          2
   *   erfc(x)  =  --------     |    exp( - t  ) dt
   *               sqrt(pi)   | |
   *                           -
   *                            x
   * </pre>
   * <b>Implementation:</b>
   * For small x, <tt>erfc(x) = 1 - erf(x)</tt>; otherwise rational
   * approximations are computed.
   * <p>
   * Code adapted from the <A HREF="http://www.sci.usq.edu.au/staff/leighb/graph/Top.html">
   * Java 2D Graph Package 2.4</A>,
   * which in turn is a port from the
   * <A HREF="http://people.ne.mediaone.net/moshier/index.html#Cephes">Cephes 2.2</A>
   * Math Library (C).
   *
   * @param a the argument to the function.
   */
  public static double errorFunctionComplemented(double a) { 
    double x,y,z,p,q;
  
    double P[] = {
      2.46196981473530512524E-10,
      5.64189564831068821977E-1,
      7.46321056442269912687E0,
      4.86371970985681366614E1,
      1.96520832956077098242E2,
      5.26445194995477358631E2,
      9.34528527171957607540E2,
      1.02755188689515710272E3,
      5.57535335369399327526E2
    };
    double Q[] = {
      //1.0
      1.32281951154744992508E1,
      8.67072140885989742329E1,
      3.54937778887819891062E2,
      9.75708501743205489753E2,
      1.82390916687909736289E3,
      2.24633760818710981792E3,
      1.65666309194161350182E3,
      5.57535340817727675546E2
    };
  
    double R[] = {
      5.64189583547755073984E-1,
      1.27536670759978104416E0,
      5.01905042251180477414E0,
      6.16021097993053585195E0,
      7.40974269950448939160E0,
      2.97886665372100240670E0
    };
    double S[] = {
      //1.00000000000000000000E0, 
      2.26052863220117276590E0,
      9.39603524938001434673E0,
      1.20489539808096656605E1,
      1.70814450747565897222E1,
      9.60896809063285878198E0,
      3.36907645100081516050E0
    };
  
    if( a < 0.0 )   x = -a;
    else            x = a;
  
    if( x < 1.0 )   return 1.0 - errorFunction(a);
  
    z = -a * a;
  
    if( z < -MAXLOG ) {
      if( a < 0 )  return( 2.0 );
      else         return( 0.0 );
    }
  
    z = Math.exp(z);
  
    if( x < 8.0 ) {
      p = polevl( x, P, 8 );
      q = p1evl( x, Q, 8 );
    } else {
      p = polevl( x, R, 5 );
      q = p1evl( x, S, 6 );
    }
  
    y = (z * p)/q;
  
    if( a < 0 ) y = 2.0 - y;
  
    if( y == 0.0 ) {
      if( a < 0 ) return 2.0;
      else        return( 0.0 );
    }
    return y;
  }
  
  /**
   * Evaluates the given polynomial of degree <tt>N</tt> at <tt>x</tt>.
   * Evaluates polynomial when coefficient of N is 1.0.
   * Otherwise same as <tt>polevl()</tt>.
   * <pre>
   *                     2          N
   * y  =  C  + C x + C x  +...+ C x
   *        0    1     2          N
   *
   * Coefficients are stored in reverse order:
   *
   * coef[0] = C  , ..., coef[N] = C  .
   *            N                   0
   * </pre>
   * The function <tt>p1evl()</tt> assumes that <tt>coef[N] = 1.0</tt> and is
   * omitted from the array.  Its calling arguments are
   * otherwise the same as <tt>polevl()</tt>.
   * <p>
   * In the interest of speed, there are no checks for out of bounds arithmetic.
   *
   * @param x argument to the polynomial.
   * @param coef the coefficients of the polynomial.
   * @param N the degree of the polynomial.
   */
  public static double p1evl( double x, double coef[], int N ) {
  
    double ans;
    ans = x + coef[0];
  
    for(int i=1; i<N; i++) ans = ans*x+coef[i];
  
    return ans;
  }

  /**
   * Evaluates the given polynomial of degree <tt>N</tt> at <tt>x</tt>.
   * <pre>
   *                     2          N
   * y  =  C  + C x + C x  +...+ C x
   *        0    1     2          N
   *
   * Coefficients are stored in reverse order:
   *
   * coef[0] = C  , ..., coef[N] = C  .
   *            N                   0
   * </pre>
   * In the interest of speed, there are no checks for out of bounds arithmetic.
   *
   * @param x argument to the polynomial.
   * @param coef the coefficients of the polynomial.
   * @param N the degree of the polynomial.
   */
  public static double polevl( double x, double coef[], int N ) {

    double ans;
    ans = coef[0];
  
    for(int i=1; i<=N; i++) ans = ans*x+coef[i];
  
    return ans;
  }

  /**
   * Returns the Incomplete Gamma function.
   * @param a the parameter of the gamma distribution.
   * @param x the integration end point.
   */
  public static double incompleteGamma(double a, double x) 
    { 
 
    double ans, ax, c, r;
 
    if( x <= 0 || a <= 0 ) return 0.0;
 
    if( x > 1.0 && x > a ) return 1.0 - incompleteGammaComplement(a,x);

    /* Compute  x**a * exp(-x) / gamma(a)  */
    ax = a * Math.log(x) - x - lnGamma(a);
    if( ax < -MAXLOG ) return( 0.0 );

    ax = Math.exp(ax);

    /* power series */
    r = a;
    c = 1.0;
    ans = 1.0;

    do {
      r += 1.0;
      c *= x/r;
      ans += c;
    }
    while( c/ans > MACHEP );
 
    return( ans * ax/a );
  }

  /**
   * Returns the Complemented Incomplete Gamma function.
   * @param a the parameter of the gamma distribution.
   * @param x the integration start point.
   */
  public static double incompleteGammaComplement( double a, double x ) {

    double ans, ax, c, yc, r, t, y, z;
    double pk, pkm1, pkm2, qk, qkm1, qkm2;

    if( x <= 0 || a <= 0 ) return 1.0;
 
    if( x < 1.0 || x < a ) return 1.0 - incompleteGamma(a,x);
 
    ax = a * Math.log(x) - x - lnGamma(a);
    if( ax < -MAXLOG ) return 0.0;
 
    ax = Math.exp(ax);
 
    /* continued fraction */
    y = 1.0 - a;
    z = x + y + 1.0;
    c = 0.0;
    pkm2 = 1.0;
    qkm2 = x;
    pkm1 = x + 1.0;
    qkm1 = z * x;
    ans = pkm1/qkm1;
 
    do {
      c += 1.0;
      y += 1.0;
      z += 2.0;
      yc = y * c;
      pk = pkm1 * z  -  pkm2 * yc;
      qk = qkm1 * z  -  qkm2 * yc;
      if( qk != 0 ) {
	r = pk/qk;
        t = Math.abs( (ans - r)/r );
	ans = r;
      } else
	t = 1.0;

      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      if( Math.abs(pk) > big ) {
	pkm2 *= biginv;
        pkm1 *= biginv;
	qkm2 *= biginv;
	qkm1 *= biginv;
      }
    } while( t > MACHEP );
 
    return ans * ax;
  }

  /**
   * Returns the Gamma function of the argument.
   */
  public static double gamma(double x) {

    double P[] = {
      1.60119522476751861407E-4,
      1.19135147006586384913E-3,
      1.04213797561761569935E-2,
      4.76367800457137231464E-2,
      2.07448227648435975150E-1,
      4.94214826801497100753E-1,
      9.99999999999999996796E-1
    };
    double Q[] = {
      -2.31581873324120129819E-5,
      5.39605580493303397842E-4,
      -4.45641913851797240494E-3,
      1.18139785222060435552E-2,
      3.58236398605498653373E-2,
      -2.34591795718243348568E-1,
      7.14304917030273074085E-2,
      1.00000000000000000320E0
    };

    double p, z;
    int i;

    double q = Math.abs(x);

    if( q > 33.0 ) {
      if( x < 0.0 ) {
	p = Math.floor(q);
	if( p == q ) throw new ArithmeticException("gamma: overflow");
	i = (int)p;
	z = q - p;
	if( z > 0.5 ) {
	  p += 1.0;
	  z = q - p;
	}
	z = q * Math.sin( Math.PI * z );
	if( z == 0.0 ) throw new ArithmeticException("gamma: overflow");
	z = Math.abs(z);
	z = Math.PI/(z * stirlingFormula(q) );

	return -z;
      } else {
	return stirlingFormula(x);
      }
    }

    z = 1.0;
    while( x >= 3.0 ) {
      x -= 1.0;
      z *= x;
    }

    while( x < 0.0 ) {
      if( x == 0.0 ) {
	throw new ArithmeticException("gamma: singular");
      } else
	if( x > -1.E-9 ) {
	  return( z/((1.0 + 0.5772156649015329 * x) * x) );
	}
      z /= x;
      x += 1.0;
    }

    while( x < 2.0 ) {
      if( x == 0.0 ) {
	throw new ArithmeticException("gamma: singular");
      } else
	if( x < 1.e-9 ) {
	  return( z/((1.0 + 0.5772156649015329 * x) * x) );
	}
      z /= x;
      x += 1.0;
    }

    if( (x == 2.0) || (x == 3.0) ) 	return z;

    x -= 2.0;
    p = polevl( x, P, 6 );
    q = polevl( x, Q, 7 );
    return  z * p / q;
  }

  /**
   * Returns the Gamma function computed by Stirling's formula.
   * The polynomial STIR is valid for 33 <= x <= 172.
   */
  public static double stirlingFormula(double x) {

    double STIR[] = {
      7.87311395793093628397E-4,
      -2.29549961613378126380E-4,
      -2.68132617805781232825E-3,
      3.47222221605458667310E-3,
      8.33333333333482257126E-2,
    };
    double MAXSTIR = 143.01608;

    double w = 1.0/x;
    double  y = Math.exp(x);

    w = 1.0 + w * polevl( w, STIR, 4 );

    if( x > MAXSTIR ) {
      /* Avoid overflow in Math.pow() */
      double v = Math.pow( x, 0.5 * x - 0.25 );
      y = v * (v / y);
    } else {
      y = Math.pow( x, x - 0.5 ) / y;
    }
    y = SQTPI * y * w;
    return y;
  }

  /**
   * Returns the Incomplete Beta Function evaluated from zero to <tt>xx</tt>.
   *
   * @param aa the alpha parameter of the beta distribution.
   * @param bb the beta parameter of the beta distribution.
   * @param xx the integration end point.
   */
  public static double incompleteBeta( double aa, double bb, double xx ) {

    double a, b, t, x, xc, w, y;
    boolean flag;

    if( aa <= 0.0 || bb <= 0.0 ) throw new 
      ArithmeticException("ibeta: Domain error!");

    if( (xx <= 0.0) || ( xx >= 1.0) ) {
      if( xx == 0.0 ) return 0.0;
      if( xx == 1.0 ) return 1.0;
      throw new ArithmeticException("ibeta: Domain error!");
    }

    flag = false;
    if( (bb * xx) <= 1.0 && xx <= 0.95) {
      t = powerSeries(aa, bb, xx);
      return t;
    }

    w = 1.0 - xx;

    /* Reverse a and b if x is greater than the mean. */
    if( xx > (aa/(aa+bb)) ) {
      flag = true;
      a = bb;
      b = aa;
      xc = xx;
      x = w;
    } else {
      a = aa;
      b = bb;
      xc = w;
      x = xx;
    }

    if( flag  && (b * x) <= 1.0 && x <= 0.95) {
      t = powerSeries(a, b, x);
      if( t <= MACHEP ) 	t = 1.0 - MACHEP;
      else  		        t = 1.0 - t;
      return t;
    }

    /* Choose expansion for better convergence. */
    y = x * (a+b-2.0) - (a-1.0);
    if( y < 0.0 )
      w = incompleteBetaFraction1( a, b, x );
    else
      w = incompleteBetaFraction2( a, b, x ) / xc;

    /* Multiply w by the factor
       a      b   _             _     _
       x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   */

    y = a * Math.log(x);
    t = b * Math.log(xc);
    if( (a+b) < MAXGAM && Math.abs(y) < MAXLOG && Math.abs(t) < MAXLOG ) {
      t = Math.pow(xc,b);
      t *= Math.pow(x,a);
      t /= a;
      t *= w;
      t *= gamma(a+b) / (gamma(a) * gamma(b));
      if( flag ) {
	if( t <= MACHEP ) 	t = 1.0 - MACHEP;
	else  		        t = 1.0 - t;
      }
      return t;
    }
    /* Resort to logarithms.  */
    y += t + lnGamma(a+b) - lnGamma(a) - lnGamma(b);
    y += Math.log(w/a);
    if( y < MINLOG )
      t = 0.0;
    else
      t = Math.exp(y);

    if( flag ) {
      if( t <= MACHEP ) 	t = 1.0 - MACHEP;
      else  		        t = 1.0 - t;
    }
    return t;
  }   

  /**
   * Continued fraction expansion #1 for incomplete beta integral.
   */
  public static double incompleteBetaFraction1( double a, double b, double x ) {

    double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
    double k1, k2, k3, k4, k5, k6, k7, k8;
    double r, t, ans, thresh;
    int n;

    k1 = a;
    k2 = a + b;
    k3 = a;
    k4 = a + 1.0;
    k5 = 1.0;
    k6 = b - 1.0;
    k7 = k4;
    k8 = a + 2.0;

    pkm2 = 0.0;
    qkm2 = 1.0;
    pkm1 = 1.0;
    qkm1 = 1.0;
    ans = 1.0;
    r = 1.0;
    n = 0;
    thresh = 3.0 * MACHEP;
    do {
      xk = -( x * k1 * k2 )/( k3 * k4 );
      pk = pkm1 +  pkm2 * xk;
      qk = qkm1 +  qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;

      xk = ( x * k5 * k6 )/( k7 * k8 );
      pk = pkm1 +  pkm2 * xk;
      qk = qkm1 +  qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;

      if( qk != 0 )		r = pk/qk;
      if( r != 0 ) {
	t = Math.abs( (ans - r)/r );
	ans = r;
      }	else
	t = 1.0;

      if( t < thresh ) return ans;

      k1 += 1.0;
      k2 += 1.0;
      k3 += 2.0;
      k4 += 2.0;
      k5 += 1.0;
      k6 -= 1.0;
      k7 += 2.0;
      k8 += 2.0;

      if( (Math.abs(qk) + Math.abs(pk)) > big ) {
	pkm2 *= biginv;
	pkm1 *= biginv;
	qkm2 *= biginv;
	qkm1 *= biginv;
      }
      if( (Math.abs(qk) < biginv) || (Math.abs(pk) < biginv) ) {
	pkm2 *= big;
	pkm1 *= big;
	qkm2 *= big;
	qkm1 *= big;
      }
    } while( ++n < 300 );

    return ans;
  }   

  /**
   * Continued fraction expansion #2 for incomplete beta integral.
   */
  public static double incompleteBetaFraction2( double a, double b, double x ) {

    double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
    double k1, k2, k3, k4, k5, k6, k7, k8;
    double r, t, ans, z, thresh;
    int n;

    k1 = a;
    k2 = b - 1.0;
    k3 = a;
    k4 = a + 1.0;
    k5 = 1.0;
    k6 = a + b;
    k7 = a + 1.0;;
    k8 = a + 2.0;

    pkm2 = 0.0;
    qkm2 = 1.0;
    pkm1 = 1.0;
    qkm1 = 1.0;
    z = x / (1.0-x);
    ans = 1.0;
    r = 1.0;
    n = 0;
    thresh = 3.0 * MACHEP;
    do {
      xk = -( z * k1 * k2 )/( k3 * k4 );
      pk = pkm1 +  pkm2 * xk;
      qk = qkm1 +  qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;

      xk = ( z * k5 * k6 )/( k7 * k8 );
      pk = pkm1 +  pkm2 * xk;
      qk = qkm1 +  qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;

      if( qk != 0 )  r = pk/qk;
      if( r != 0 ) {
	t = Math.abs( (ans - r)/r );
	ans = r;
      } else
	t = 1.0;

      if( t < thresh ) return ans;

      k1 += 1.0;
      k2 -= 1.0;
      k3 += 2.0;
      k4 += 2.0;
      k5 += 1.0;
      k6 += 1.0;
      k7 += 2.0;
      k8 += 2.0;

      if( (Math.abs(qk) + Math.abs(pk)) > big ) {
	pkm2 *= biginv;
	pkm1 *= biginv;
	qkm2 *= biginv;
	qkm1 *= biginv;
      }
      if( (Math.abs(qk) < biginv) || (Math.abs(pk) < biginv) ) {
	pkm2 *= big;
	pkm1 *= big;
	qkm2 *= big;
	qkm1 *= big;
      }
    } while( ++n < 300 );

    return ans;
  }

  /**
   * Power series for incomplete beta integral.
   * Use when b*x is small and x not too close to 1.  
   */
  public static double powerSeries( double a, double b, double x ) {

    double s, t, u, v, n, t1, z, ai;
    
    ai = 1.0 / a;
    u = (1.0 - b) * x;
    v = u / (a + 1.0);
    t1 = v;
    t = u;
    n = 2.0;
    s = 0.0;
    z = MACHEP * ai;
    while( Math.abs(v) > z ) {
      u = (n - b) * x / n;
      t *= u;
      v = t / (a + n);
      s += v; 
      n += 1.0;
    }
    s += t1;
    s += ai;

    u = a * Math.log(x);
    if( (a+b) < MAXGAM && Math.abs(u) < MAXLOG ) {
      t = gamma(a+b)/(gamma(a)*gamma(b));
      s = s * t * Math.pow(x,a);
    } else {
      t = lnGamma(a+b) - lnGamma(a) - lnGamma(b) + u + Math.log(s);
      if( t < MINLOG ) 	s = 0.0;
      else  	            s = Math.exp(t);
    }
    return s;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5616 $");
  }

  /**
   * Main method for testing this class.
   */
  public static void main(String[] ops) {

    System.out.println("Binomial standard error (0.5, 100): " + 
		       Statistics.binomialStandardError(0.5, 100));
    System.out.println("Chi-squared probability (2.558, 10): " +
		       Statistics.chiSquaredProbability(2.558, 10));
    System.out.println("Normal probability (0.2): " +
		       Statistics.normalProbability(0.2));
    System.out.println("F probability (5.1922, 4, 5): " +
		       Statistics.FProbability(5.1922, 4, 5));
    System.out.println("lnGamma(6): "+ Statistics.lnGamma(6));
  }
}
