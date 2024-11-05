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
 *    PoissonEstimator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;

import weka.core.Capabilities.Capability;
import weka.core.Capabilities;
import weka.core.RevisionUtils;
import weka.core.Utils;

/** 
 * Simple probability estimator that places a single Poisson distribution
 * over the observed values.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 5490 $
 */
public class PoissonEstimator
  extends Estimator
  implements IncrementalEstimator {

  /** for serialization */
  private static final long serialVersionUID = 7669362595289236662L;
  
  /** The number of values seen */
  private double m_NumValues;
  
  /** The sum of the values seen */
  private double m_SumOfValues;
  
  /** 
   * The average number of times
   * an event occurs in an interval.
   */
  private double m_Lambda;
  
  
  /**
   * Calculates the log factorial of a number.
   *
   * @param x input number.
   * @return log factorial of x.
   */
  private double logFac(double x) {
    
    double result = 0;
    for (double i = 2; i <= x; i++) {
      result += Math.log(i);
    }
    return result;
  }
  
  /**
   * Returns value for Poisson distribution
   *
   * @param x the argument to the kernel function
   * @return the value for a Poisson kernel
   */
  private double Poisson(double x) {
    
    return Math.exp(-m_Lambda + (x * Math.log(m_Lambda)) - logFac(x));
  }
  
  /**
   * Add a new data value to the current estimator.
   *
   * @param data the new data value 
   * @param weight the weight assigned to the data value 
   */
  public void addValue(double data, double weight) {
    
    m_NumValues += weight;
    m_SumOfValues += data * weight;
    if (m_NumValues != 0) {
      m_Lambda = m_SumOfValues / m_NumValues;
    }
  }
  
  /**
   * Get a probability estimate for a value
   *
   * @param data the value to estimate the probability of
   * @return the estimated probability of the supplied value
   */
  public double getProbability(double data) {
    
    return Poisson(data);
  }
  
  /** Display a representation of this estimator */
  public String toString() {
    
    return "Poisson Lambda = " + Utils.doubleToString(m_Lambda, 4, 2) + "\n";
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
      if (argv.length == 0) {
        System.out.println("Please specify a set of instances.");
        return;
      }
      PoissonEstimator newEst = new PoissonEstimator();
      for(int i = 0; i < argv.length; i++) {
        double current = Double.valueOf(argv[i]).doubleValue();
        System.out.println(newEst);
        System.out.println("Prediction for " + current 
            + " = " + newEst.getProbability(current));
        newEst.addValue(current, 1);
      }
      
    } catch (Exception e) {
      System.out.println(e.getMessage());
    }
  }
}
