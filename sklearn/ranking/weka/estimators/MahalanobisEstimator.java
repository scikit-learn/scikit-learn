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
 *    MahalanobisEstimator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;

import weka.core.Capabilities.Capability;
import weka.core.matrix.Matrix;
import weka.core.Capabilities;
import weka.core.RevisionUtils;
import weka.core.Utils;

/** 
 * Simple probability estimator that places a single normal distribution
 * over the observed values.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 5490 $
 */
public class MahalanobisEstimator extends Estimator implements IncrementalEstimator {
  
  /** for serialization  */
  private static final long serialVersionUID = 8950225468990043868L;
  
  /** The inverse of the covariance matrix */
  private Matrix m_CovarianceInverse;
  
  /** The determinant of the covariance matrix */
  private double m_Determinant;
  
  /**
   * The difference between the conditioning value and the conditioning mean
   */
  private double m_ConstDelta;
  
  /** The mean of the values */
  private double m_ValueMean;
  
  /** 2 * PI */
  private static double TWO_PI = 2 * Math.PI;
  
  /**
   * Returns value for normal kernel
   *
   * @param x the argument to the kernel function
   * @param variance the variance
   * @return the value for a normal kernel
   */
  private double normalKernel(double x) {
    
    Matrix thisPoint = new Matrix(1, 2);
    thisPoint.set(0, 0, x);
    thisPoint.set(0, 1, m_ConstDelta);
    return Math.exp(-thisPoint.times(m_CovarianceInverse).
        times(thisPoint.transpose()).get(0, 0) 
        / 2) / (Math.sqrt(TWO_PI) * m_Determinant);
  }
  
  /**
   * Constructor
   *
   * @param covariance
   * @param constDelta
   * @param valueMean
   */
  public MahalanobisEstimator(Matrix covariance, double constDelta,
      double valueMean) {
    
    m_CovarianceInverse = null;
    if ((covariance.getRowDimension() == 2) && (covariance.getColumnDimension() == 2)) {
      double a = covariance.get(0, 0);
      double b = covariance.get(0, 1);
      double c = covariance.get(1, 0);
      double d = covariance.get(1, 1);
      if (a == 0) {
        a = c; c = 0;
        double temp = b;
        b = d; d = temp;
      }
      if (a == 0) {
        return;
      }
      double denom = d - c * b / a;
      if (denom == 0) {
        return;
      }
      m_Determinant = covariance.get(0, 0) * covariance.get(1, 1)
      - covariance.get(1, 0) * covariance.get(0, 1);
      m_CovarianceInverse = new Matrix(2, 2);
      m_CovarianceInverse.set(0, 0, 1.0 / a + b * c / a / a / denom);
      m_CovarianceInverse.set(0, 1, -b / a / denom);
      m_CovarianceInverse.set(1, 0, -c / a / denom);
      m_CovarianceInverse.set(1, 1, 1.0 / denom);
      m_ConstDelta = constDelta;
      m_ValueMean = valueMean;
    }
  }
  
  /**
   * Add a new data value to the current estimator. Does nothing because the
   * data is provided in the constructor.
   *
   * @param data the new data value 
   * @param weight the weight assigned to the data value 
   */
  public void addValue(double data, double weight) {
    
  }
  
  /**
   * Get a probability estimate for a value
   *
   * @param data the value to estimate the probability of
   * @return the estimated probability of the supplied value
   */
  public double getProbability(double data) {
    
    double delta = data - m_ValueMean;
    if (m_CovarianceInverse == null) {
      return 0;
    }
    return normalKernel(delta);
  }
  
  /** Display a representation of this estimator */
  public String toString() {
    
    if (m_CovarianceInverse == null) {
      return "No covariance inverse\n";
    }
    return "Mahalanovis Distribution. Mean = "
    + Utils.doubleToString(m_ValueMean, 4, 2)
    + "  ConditionalOffset = "
    + Utils.doubleToString(m_ConstDelta, 4, 2) + "\n"
    + "Covariance Matrix: Determinant = " + m_Determinant 
    + "  Inverse:\n" + m_CovarianceInverse;
  }
  
  /**
   * Returns default capabilities of the classifier.
   *
   * @return      the capabilities of this classifier
   */
  public Capabilities getCapabilities() {
    Capabilities result = super.getCapabilities();
    result.disableAll();
    
    // class
    if (!m_noClass) {
      result.enable(Capability.NOMINAL_CLASS);
      result.enable(Capability.MISSING_CLASS_VALUES);
    } else {
      result.enable(Capability.NO_CLASS);
    }
    
    // attributes
    result.enable(Capability.NUMERIC_ATTRIBUTES);
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 5490 $");
  }
  
  /**
   * Main method for testing this class.
   *
   * @param argv should contain a sequence of numeric values
   */
  public static void main(String [] argv) {
    
    try {
      double delta = 0.5;
      double xmean = 0;
      double lower = 0;
      double upper = 10;
      Matrix covariance = new Matrix(2, 2);
      covariance.set(0, 0, 2);
      covariance.set(0, 1, -3);
      covariance.set(1, 0, -4);
      covariance.set(1, 1, 5);
      if (argv.length > 0) {
        covariance.set(0, 0, Double.valueOf(argv[0]).doubleValue());
      }
      if (argv.length > 1) {
        covariance.set(0, 1, Double.valueOf(argv[1]).doubleValue());
      }
      if (argv.length > 2) {
        covariance.set(1, 0, Double.valueOf(argv[2]).doubleValue());
      }
      if (argv.length > 3) {
        covariance.set(1, 1, Double.valueOf(argv[3]).doubleValue());
      }
      if (argv.length > 4) {
        delta = Double.valueOf(argv[4]).doubleValue();
      }
      if (argv.length > 5) {
        xmean = Double.valueOf(argv[5]).doubleValue();
      }
      
      MahalanobisEstimator newEst = new MahalanobisEstimator(covariance,
          delta, xmean);
      if (argv.length > 6) {
        lower = Double.valueOf(argv[6]).doubleValue();
        if (argv.length > 7) {
          upper = Double.valueOf(argv[7]).doubleValue();
        }
        double increment = (upper - lower) / 50;
        for(double current = lower; current <= upper; current+= increment)
          System.out.println(current + "  " + newEst.getProbability(current));
      } else {
        System.out.println("Covariance Matrix\n" + covariance);
        System.out.println(newEst);
      }
    } catch (Exception e) {
      System.out.println(e.getMessage());
    }
  }
}
