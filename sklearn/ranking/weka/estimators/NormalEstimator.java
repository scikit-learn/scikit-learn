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
 *    NormalEstimator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;

import weka.core.Capabilities.Capability;
import weka.core.Capabilities;
import weka.core.RevisionUtils;
import weka.core.Statistics;
import weka.core.Utils;

/** 
 * Simple probability estimator that places a single normal distribution
 * over the observed values.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 5490 $
 */
public class NormalEstimator
  extends Estimator
  implements IncrementalEstimator {

  /** for serialization */
  private static final long serialVersionUID = 93584379632315841L;

  /** The sum of the weights */
  private double m_SumOfWeights;

  /** The sum of the values seen */
  private double m_SumOfValues;

  /** The sum of the values squared */
  private double m_SumOfValuesSq;

  /** The current mean */
  private double m_Mean;

  /** The current standard deviation */
  private double m_StandardDev;

  /** The precision of numeric values ( = minimum std dev permitted) */
  private double m_Precision;

  /**
   * Round a data value using the defined precision for this estimator
   *
   * @param data the value to round
   * @return the rounded data value
   */
  private double round(double data) {

    return Math.rint(data / m_Precision) * m_Precision;
  }
  
  // ===============
  // Public methods.
  // ===============
  
  /**
   * Constructor that takes a precision argument.
   *
   * @param precision the precision to which numeric values are given. For
   * example, if the precision is stated to be 0.1, the values in the
   * interval (0.25,0.35] are all treated as 0.3. 
   */
  public NormalEstimator(double precision) {

    m_Precision = precision;

    // Allow at most 3 sd's within one interval
    m_StandardDev = m_Precision / (2 * 3);
  }

  /**
   * Add a new data value to the current estimator.
   *
   * @param data the new data value 
   * @param weight the weight assigned to the data value 
   */
  public void addValue(double data, double weight) {

    if (weight == 0) {
      return;
    }
    data = round(data);
    m_SumOfWeights += weight;
    m_SumOfValues += data * weight;
    m_SumOfValuesSq += data * data * weight;

    if (m_SumOfWeights > 0) {
      m_Mean = m_SumOfValues / m_SumOfWeights;
      double stdDev = Math.sqrt(Math.abs(m_SumOfValuesSq 
					  - m_Mean * m_SumOfValues) 
					 / m_SumOfWeights);
      // If the stdDev ~= 0, we really have no idea of scale yet, 
      // so stick with the default. Otherwise...
      if (stdDev > 1e-10) {
	m_StandardDev = Math.max(m_Precision / (2 * 3), 
				 // allow at most 3sd's within one interval 
				 stdDev);
      }
    }
  }

  /**
   * Get a probability estimate for a value
   *
   * @param data the value to estimate the probability of
   * @return the estimated probability of the supplied value
   */
  public double getProbability(double data) {

    data = round(data);
    double zLower = (data - m_Mean - (m_Precision / 2)) / m_StandardDev;
    double zUpper = (data - m_Mean + (m_Precision / 2)) / m_StandardDev;
    
    double pLower = Statistics.normalProbability(zLower);
    double pUpper = Statistics.normalProbability(zUpper);
    return pUpper - pLower;
  }

  /**
   * Display a representation of this estimator
   */
  public String toString() {

    return "Normal Distribution. Mean = " + Utils.doubleToString(m_Mean, 4)
      + " StandardDev = " + Utils.doubleToString(m_StandardDev, 4)
      + " WeightSum = " + Utils.doubleToString(m_SumOfWeights, 4)
      + " Precision = " + m_Precision + "\n";
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
   * Return the value of the mean of this normal estimator.
   *
   * @return the mean
   */
  public double getMean() {
    return m_Mean;
  }

  /**
   * Return the value of the standard deviation of this normal estimator.
   *
   * @return the standard deviation
   */
  public double getStdDev() {
    return m_StandardDev;
  }

  /**
   * Return the value of the precision of this normal estimator.
   *
   * @return the precision
   */
  public double getPrecision() {
    return m_Precision;
  }

  /**
   * Return the sum of the weights for this normal estimator.
   *
   * @return the sum of the weights
   */
  public double getSumOfWeights() {
    return m_SumOfWeights;
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
      NormalEstimator newEst = new NormalEstimator(0.01);
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
