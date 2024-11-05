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
 *    KernelEstimator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;

import weka.core.Capabilities.Capability;
import weka.core.Capabilities;
import weka.core.RevisionUtils;
import weka.core.Utils;
import weka.core.Statistics;

/** 
 * Simple kernel density estimator. Uses one gaussian kernel per observed
 * data value.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 5490 $
 */
public class KernelEstimator extends Estimator implements IncrementalEstimator {

  /** for serialization */
  private static final long serialVersionUID = 3646923563367683925L;

  /** Vector containing all of the values seen */
  private double [] m_Values;

  /** Vector containing the associated weights */
  private double [] m_Weights;

  /** Number of values stored in m_Weights and m_Values so far */
  private int m_NumValues;

  /** The sum of the weights so far */
  private double m_SumOfWeights;

  /** The standard deviation */
  private double m_StandardDev;

  /** The precision of data values */
  private double m_Precision;

  /** Whether we can optimise the kernel summation */
  private boolean m_AllWeightsOne;

  /** Maximum percentage error permitted in probability calculations */
  private static double MAX_ERROR = 0.01;


  /**
   * Execute a binary search to locate the nearest data value
   *
   * @param the data value to locate
   * @return the index of the nearest data value
   */
  private int findNearestValue(double key) {

    int low = 0; 
    int high = m_NumValues;
    int middle = 0;
    while (low < high) {
      middle = (low + high) / 2;
      double current = m_Values[middle];
      if (current == key) {
	return middle;
      }
      if (current > key) {
	high = middle;
      } else if (current < key) {
	low = middle + 1;
      }
    }
    return low;
  }

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
   * @param precision the  precision to which numeric values are given. For
   * example, if the precision is stated to be 0.1, the values in the
   * interval (0.25,0.35] are all treated as 0.3. 
   */
  public KernelEstimator(double precision) {

    m_Values = new double [50];
    m_Weights = new double [50];
    m_NumValues = 0;
    m_SumOfWeights = 0;
    m_AllWeightsOne = true;
    m_Precision = precision;
    // precision cannot be zero
    if (m_Precision < Utils.SMALL) m_Precision = Utils.SMALL;
    //    m_StandardDev = 1e10 * m_Precision; // Set the standard deviation initially very wide
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
    int insertIndex = findNearestValue(data);
    if ((m_NumValues <= insertIndex) || (m_Values[insertIndex] != data)) {
      if (m_NumValues < m_Values.length) {
        int left = m_NumValues - insertIndex; 
        System.arraycopy(m_Values, insertIndex, 
            m_Values, insertIndex + 1, left);
        System.arraycopy(m_Weights, insertIndex, 
            m_Weights, insertIndex + 1, left);
        
        m_Values[insertIndex] = data;
        m_Weights[insertIndex] = weight;
        m_NumValues++;
      } else {
        double [] newValues = new double [m_Values.length * 2];
        double [] newWeights = new double [m_Values.length * 2];
        int left = m_NumValues - insertIndex; 
        System.arraycopy(m_Values, 0, newValues, 0, insertIndex);
        System.arraycopy(m_Weights, 0, newWeights, 0, insertIndex);
        newValues[insertIndex] = data;
        newWeights[insertIndex] = weight;
        System.arraycopy(m_Values, insertIndex, 
            newValues, insertIndex + 1, left);
        System.arraycopy(m_Weights, insertIndex, 
            newWeights, insertIndex + 1, left);
        m_NumValues++;
        m_Values = newValues;
        m_Weights = newWeights;
      }
      if (weight != 1) {
        m_AllWeightsOne = false;
      }
    } else {
      m_Weights[insertIndex] += weight;
      m_AllWeightsOne = false;      
    }
    m_SumOfWeights += weight;
    double range = m_Values[m_NumValues - 1] - m_Values[0];
    if (range > 0) {
      m_StandardDev = Math.max(range / Math.sqrt(m_SumOfWeights), 
          // allow at most 3 sds within one interval
          m_Precision / (2 * 3));
    }
  }
  
  /**
   * Get a probability estimate for a value.
   *
   * @param data the value to estimate the probability of
   * @return the estimated probability of the supplied value
   */
  public double getProbability(double data) {

    double delta = 0, sum = 0, currentProb = 0;
    double zLower = 0, zUpper = 0;
    if (m_NumValues == 0) {
      zLower = (data - (m_Precision / 2)) / m_StandardDev;
      zUpper = (data + (m_Precision / 2)) / m_StandardDev;
      return (Statistics.normalProbability(zUpper)
	      - Statistics.normalProbability(zLower));
    }
    double weightSum = 0;
    int start = findNearestValue(data);
    for (int i = start; i < m_NumValues; i++) {
      delta = m_Values[i] - data;
      zLower = (delta - (m_Precision / 2)) / m_StandardDev;
      zUpper = (delta + (m_Precision / 2)) / m_StandardDev;
      currentProb = (Statistics.normalProbability(zUpper)
		     - Statistics.normalProbability(zLower));
      sum += currentProb * m_Weights[i];
      /*
      System.out.print("zL" + (i + 1) + ": " + zLower + " ");
      System.out.print("zU" + (i + 1) + ": " + zUpper + " ");
      System.out.print("P" + (i + 1) + ": " + currentProb + " ");
      System.out.println("total: " + (currentProb * m_Weights[i]) + " ");
      */
      weightSum += m_Weights[i];
      if (currentProb * (m_SumOfWeights - weightSum) < sum * MAX_ERROR) {
	break;
      }
    }
    for (int i = start - 1; i >= 0; i--) {
      delta = m_Values[i] - data;
      zLower = (delta - (m_Precision / 2)) / m_StandardDev;
      zUpper = (delta + (m_Precision / 2)) / m_StandardDev;
      currentProb = (Statistics.normalProbability(zUpper)
		     - Statistics.normalProbability(zLower));
      sum += currentProb * m_Weights[i];
      weightSum += m_Weights[i];
      if (currentProb * (m_SumOfWeights - weightSum) < sum * MAX_ERROR) {
	break;
      }
    }
    return sum / m_SumOfWeights;
  }

  /** Display a representation of this estimator */
  public String toString() {

    String result = m_NumValues + " Normal Kernels. \nStandardDev = " 
      + Utils.doubleToString(m_StandardDev,6,4)
      + " Precision = " + m_Precision;
    if (m_NumValues == 0) {
      result += "  \nMean = 0";
    } else {
      result += "  \nMeans =";
      for (int i = 0; i < m_NumValues; i++) {
	result += " " + m_Values[i];
      }
      if (!m_AllWeightsOne) {
	result += "\nWeights = ";
	for (int i = 0; i < m_NumValues; i++) {
	  result += " " + m_Weights[i];
	}
      }
    }
    return result + "\n";
  }

  /**
   * Return the number of kernels in this kernel estimator
   *
   * @return the number of kernels
   */
  public int getNumKernels() {
    return m_NumValues;
  }

  /**
   * Return the means of the kernels.
   *
   * @return the means of the kernels
   */
  public double[] getMeans() {
    return m_Values;
  }

  /**
   * Return the weights of the kernels.
   *
   * @return the weights of the kernels
   */
  public double[] getWeights() {
    return m_Weights;
  }

  /**
   * Return the precision of this kernel estimator.
   *
   * @return the precision
   */
  public double getPrecision() {
    return m_Precision;
  }

  /**
   * Return the standard deviation of this kernel estimator.
   *
   * @return the standard deviation
   */
  public double getStdDev() {
    return m_StandardDev;
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
      if (argv.length < 2) {
	System.out.println("Please specify a set of instances.");
	return;
      }
      KernelEstimator newEst = new KernelEstimator(0.01);
      for (int i = 0; i < argv.length - 3; i += 2) {
	newEst.addValue(Double.valueOf(argv[i]).doubleValue(), 
			Double.valueOf(argv[i + 1]).doubleValue());
      }
      System.out.println(newEst);

      double start = Double.valueOf(argv[argv.length - 2]).doubleValue();
      double finish = Double.valueOf(argv[argv.length - 1]).doubleValue();
      for (double current = start; current < finish; 
	  current += (finish - start) / 50) {
	System.out.println("Data: " + current + " " 
			   + newEst.getProbability(current));
      }
    } catch (Exception e) {
      System.out.println(e.getMessage());
    }
  }
}
