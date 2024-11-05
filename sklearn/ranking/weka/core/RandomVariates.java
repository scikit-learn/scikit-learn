/*
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 *    RandomVariates.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.util.Random;

/**
 * Class implementing some simple random variates generator.
 *
 * @author Xin Xu (xx5@cs.waikato.ac.nz)
 * @version $Revision: 5359 $
 */
public final class RandomVariates extends Random implements RevisionHandler {

  /** for serialization */
  private static final long serialVersionUID = -4763742718209460354L;

  /** 
   * Simply the constructor of super class
   */
  public RandomVariates(){ super(); }

  /** 
   * Simply the constructor of super class
   *
   * @param seed the seed in this random object
   */
  public RandomVariates(long seed){ super(seed); }

  /** 
   * Simply use the method of the super class
   *
   * @param bits - random bits
   * @return the next pseudorandom value from this random number 
   * generator's sequence.
   */
  protected int next(int bits) {return super.next(bits);}

  /**
   * Generate a value of a variate following standard exponential
   * distribution using simple inverse method.<p>
   *
   * Variates related to standard Exponential can be generated using simple
   * transformations.
   * @return a value of the variate
   */ 
  public double nextExponential(){
    return -Math.log(1.0-super.nextDouble());
  }

  /**
   * Generate a value of a variate following standard Erlang distribution.
   * It can be used when the shape parameter is an integer and not too large,
   * say, <100.  When the parameter is not an integer (which is often named
   * Gamma distribution) or is large, use <code>nextGamma(double a)</code>
   * instead.
   *
   * @param a the shape parameter, must be no less than 1
   * @return a value of the variate
   * @exception Exception if parameter less than 1
   */
  public double nextErlang(int a) throws Exception{
    if(a<1)
      throw new Exception("Shape parameter of Erlang distribution must be greater than 1!");

    double product = 1.0;
    for(int i=1; i<=a; i++)
      product *= super.nextDouble();

    return -Math.log(product);
  }

  /**
   * Generate a value of a variate following standard Gamma distribution 
   * with shape parameter a.<p>
   * If a>1, it uses a rejection method developed by Minh(1988)"Generating
   * Gamma Variates", ACM Trans. on Math. Software, Vol.14, No.3, pp261-266.
   * <br>
   * If a<1, it uses the algorithm "GS" developed by Ahrens and Dieter(1974)
   * "COMPUTER METHODS FOR SAMPLING FROM GAMMA, BETA, POISSON AND BINOMIAL
   * DISTRIBUTIONS", COMPUTING, 12 (1974), pp223-246, and further implemented
   * in Fortran by Ahrens, Kohrt and Dieter(1983) "Algorithm 599: sampling
   * from Gamma and Poisson distributions", ACM Trans. on Math. Software, 
   * Vol.9 No.2, pp255-257.<p> 
   * 
   * Variates related to standard Gamma can be generated using simple
   * transformations.
   *
   * @param a the shape parameter, must be greater than 1
   * @return a value of the variate
   * @exception Exception if parameter not greater than 1
   */
  public double nextGamma(double a) throws Exception{
    if(a<=0.0)
      throw new Exception("Shape parameter of Gamma distribution"+
      "must be greater than 0!");
    else if (a==1.0)
      return nextExponential();
    else if (a<1.0){
      double b=1.0+Math.exp(-1.0)*a, p,x, condition;
      do{
        p=b*super.nextDouble();
        if(p<1.0){
          x = Math.exp(Math.log(p)/a);
          condition = x;
        }
        else{
          x = -Math.log((b-p)/a);
          condition = (1.0-a)*Math.log(x);
        }
      }
      while(nextExponential() < condition);
      return x;	    
    }
    else{ // a>1
      double b=a-1.0, D=Math.sqrt(b), D1,x1,x2,xl,f1,f2,x4,x5,xr,f4,f5,
      p1,p2,p3,p4;

      // Initialization
      if(a<=2.0){
        D1 = b/2.0;
        x1 = 0.0;
        x2 = D1;
        xl = -1.0;
        f1 = 0.0;
      }
      else{
        D1 = D-0.5;
        x2 = b-D1;
        x1 = x2-D1;
        xl = 1.0-b/x1;
        f1 = Math.exp(b*Math.log(x1/b)+2.0*D1);
      }

      f2=Math.exp(b*Math.log(x2/b)+D1);
      x4 = b+D;
      x5 = x4+D;
      xr = 1.0-b/x5;
      f4 = Math.exp(b*Math.log(x4/b)-D);
      f5 = Math.exp(b*Math.log(x5/b)-2.0*D);
      p1 = 2.0*f4*D;
      p2 = 2.0*f2*D1+p1;
      p3 = f5/xr+p2;
      p4 = -f1/xl+p3;

      // Generation
      double u, w=Double.MAX_VALUE, x=b, v, xp;
      while(Math.log(w) > (b*Math.log(x/b)+b-x) || x < 0.0){
        u=super.nextDouble()*p4;
        if(u<=p1){ // step 5-6
          w = u/D-f4;
          if(w<=0.0) return (b+u/f4);
          if(w<=f5)  return (x4+(w*D)/f5);

          v = super.nextDouble();
          x=x4+v*D;
          xp=2.0*x4-x;

          if(w >= f4+(f4-1)*(x-x4)/(x4-b))
            return xp;
          if(w <= f4+(b/x4-1)*f4*(x-x4))
            return x;
          if((w < 2.0*f4-1.0) || 
              (w < 2.0*f4-Math.exp(b*Math.log(xp/b)+b-xp)))
            continue;
          return xp;
        }
        else if(u<=p2){ // step 7-8
          w = (u-p1)/D1-f2;
          if(w<=0.0) return (b-(u-p1)/f2);
          if(w<=f1)  return (x1+w*D1/f1);

          v = super.nextDouble();
          x=x1+v*D1;
          xp=2.0*x2-x;

          if(w >= f2+(f2-1)*(x-x2)/(x2-b))
            return xp;
          if(w <= f2*(x-x1)/D1)
            return x;
          if((w < 2.0*f2-1.0) || 
              (w < 2.0*f2-Math.exp(b*Math.log(xp/b)+b-xp)))
            continue;
          return xp;
        }
        else if(u<p3){ // step 9
          w = super.nextDouble();
          u = (p3-u)/(p3-p2);
          x = x5-Math.log(u)/xr;
          if(w <= (xr*(x5-x)+1.0)/u) return x;
          w = w*f5*u;
        }
        else{ // step 10
          w = super.nextDouble();
          u = (p4-u)/(p4-p3);
          x = x1-Math.log(u)/xl;
          if(x<0.0) continue;
          if(w <= (xl*(x1-x)+1.0)/u) return x;
          w = w*f1*u;
        }
      }

      return x;
    }	
  }

  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5359 $");
  }


  /**
   * Main method for testing this class.
   *
   * @param ops # of variates/seed, default is 10/45
   */
  public static void main(String[] ops) {

    int n = Integer.parseInt(ops[0]);
    if(n<=0)
      n=10;
    long seed = Long.parseLong(ops[1]);
    if(seed <= 0)
      seed = 45;
    RandomVariates var = new RandomVariates(seed);
    double varb[] = new double[n];

    try {
      System.out.println("Generate "+n+" values with std. exp dist:");
      for(int i=0; i<n; i++){
        varb[i] = var.nextExponential();
        System.out.print("["+i+"] "+varb[i]+", ");
      }

      System.out.println("\nMean is "+Utils.mean(varb)+
          ", Variance is "+Utils.variance(varb)+
          "\n\nGenerate "+n+" values with"+
      " std. Erlang-5 dist:");

      for(int i=0; i<n; i++){
        varb[i] = var.nextErlang(5);
        System.out.print("["+i+"] "+varb[i]+", ");
      }

      System.out.println("\nMean is "+Utils.mean(varb)+
          ", Variance is "+Utils.variance(varb)+
          "\n\nGenerate "+n+" values with"+
      " std. Gamma(4.5) dist:");

      for(int i=0; i<n; i++){
        varb[i] = var.nextGamma(4.5);
        System.out.print("["+i+"] "+varb[i]+", ");
      }	 

      System.out.println("\nMean is "+Utils.mean(varb)+
          ", Variance is "+Utils.variance(varb)+
          "\n\nGenerate "+n+" values with"+
      " std. Gamma(0.5) dist:");

      for(int i=0; i<n; i++){
        varb[i] = var.nextGamma(0.5);
        System.out.print("["+i+"] "+varb[i]+", ");
      }	  	  

      System.out.println("\nMean is "+Utils.mean(varb)+
          ", Variance is "+Utils.variance(varb)+
          "\n\nGenerate "+n+" values with"+
      " std. Gaussian(5, 2) dist:");

      for(int i=0; i<n; i++){
        varb[i] = var.nextGaussian()*2.0+5.0;
        System.out.print("["+i+"] "+varb[i]+", ");
      }	  	  
      System.out.println("\nMean is "+Utils.mean(varb)+
          ", Variance is "+Utils.variance(varb)+"\n");

    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
