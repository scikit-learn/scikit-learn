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
 *    KKConditionalEstimator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;

import java.util.Random;

import weka.core.RevisionUtils;
import weka.core.Statistics;
import weka.core.Utils;

/** 
 * Conditional probability estimator for a numeric domain conditional upon
 * a numeric domain.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.8 $
 */
public class KKConditionalEstimator implements ConditionalEstimator {

  /** Vector containing all of the values seen */
  private double [] m_Values;

  /** Vector containing all of the conditioning values seen */
  private double [] m_CondValues;

  /** Vector containing the associated weights */
  private double [] m_Weights;

  /**
   * Number of values stored in m_Weights, m_CondValues, and m_Values so far
   */
  private int m_NumValues;

  /** The sum of the weights so far */
  private double m_SumOfWeights;

  /** Current standard dev */
  private double m_StandardDev;

  /** Whether we can optimise the kernel summation */
  private boolean m_AllWeightsOne;

  /** The numeric precision */
  private double m_Precision;

  /**
   * Execute a binary search to locate the nearest data value
   *
   * @param key the data value to locate
   * @param secondaryKey the data value to locate
   * @return the index of the nearest data value
   */
  private int findNearestPair(double key, double secondaryKey) {

    int low = 0; 
    int high = m_NumValues;
    int middle = 0;
    while (low < high) {
      middle = (low + high) / 2;
      double current = m_CondValues[middle];
      if (current == key) {
	double secondary = m_Values[middle];
	if (secondary == secondaryKey) {
	  return middle;
	}
	if (secondary > secondaryKey) {
	  high = middle;
	} else if (secondary < secondaryKey) {
	  low = middle+1;
	}
      }
      if (current > key) {
	high = middle;
      } else if (current < key) {
	low = middle+1;
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
  
  /**
   * Constructor
   *
   * @param precision the  precision to which numeric values are given. For
   * example, if the precision is stated to be 0.1, the values in the
   * interval (0.25,0.35] are all treated as 0.3. 
   */
  public KKConditionalEstimator(double precision) {

    m_CondValues = new double [50];
    m_Values = new double [50];
    m_Weights = new double [50];
    m_NumValues = 0;
    m_SumOfWeights = 0;
    m_StandardDev = 0;
    m_AllWeightsOne = true;
    m_Precision = precision;
  }

  /**
   * Add a new data value to the current estimator.
   *
   * @param data the new data value 
   * @param given the new value that data is conditional upon 
   * @param weight the weight assigned to the data value 
   */
  public void addValue(double data, double given, double weight) {

    data = round(data);
    given = round(given);
    int insertIndex = findNearestPair(given, data);
    if ((m_NumValues <= insertIndex)
	|| (m_CondValues[insertIndex] != given)
	|| (m_Values[insertIndex] != data)) {
      if (m_NumValues < m_Values.length) {
	int left = m_NumValues - insertIndex; 
	System.arraycopy(m_Values, insertIndex, 
			 m_Values, insertIndex + 1, left);
	System.arraycopy(m_CondValues, insertIndex, 
			 m_CondValues, insertIndex + 1, left);
	System.arraycopy(m_Weights, insertIndex, 
			 m_Weights, insertIndex + 1, left);
	m_Values[insertIndex] = data;
	m_CondValues[insertIndex] = given;
	m_Weights[insertIndex] = weight;
	m_NumValues++;
      } else {
	double [] newValues = new double [m_Values.length*2];
	double [] newCondValues = new double [m_Values.length*2];
	double [] newWeights = new double [m_Values.length*2];
	int left = m_NumValues - insertIndex; 
	System.arraycopy(m_Values, 0, newValues, 0, insertIndex);
	System.arraycopy(m_CondValues, 0, newCondValues, 0, insertIndex);
	System.arraycopy(m_Weights, 0, newWeights, 0, insertIndex);
	newValues[insertIndex] = data;
	newCondValues[insertIndex] = given;
	newWeights[insertIndex] = weight;
	System.arraycopy(m_Values, insertIndex, 
			 newValues, insertIndex+1, left);
	System.arraycopy(m_CondValues, insertIndex, 
			 newCondValues, insertIndex+1, left);
	System.arraycopy(m_Weights, insertIndex, 
			 newWeights, insertIndex+1, left);
	m_NumValues++;
	m_Values = newValues;
	m_CondValues = newCondValues;
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
    double range = m_CondValues[m_NumValues-1] - m_CondValues[0];
    m_StandardDev = Math.max(range / Math.sqrt(m_SumOfWeights), 
			     // allow at most 3 sds within one interval
			     m_Precision / (2 * 3));
  }

  /**
   * Get a probability estimator for a value
   *
   * @param given the new value that data is conditional upon 
   * @return the estimator for the supplied value given the condition
   */
  public Estimator getEstimator(double given) {

    Estimator result = new KernelEstimator(m_Precision);
    if (m_NumValues == 0) {
      return result;
    }

    double delta = 0, currentProb = 0;
    double zLower, zUpper;
    for(int i = 0; i < m_NumValues; i++) {
      delta = m_CondValues[i] - given;
      zLower = (delta - (m_Precision / 2)) / m_StandardDev;
      zUpper = (delta + (m_Precision / 2)) / m_StandardDev;
      currentProb = (Statistics.normalProbability(zUpper)
		     - Statistics.normalProbability(zLower));
      result.addValue(m_Values[i], currentProb * m_Weights[i]);
    }
    return result;
  }

  /**
   * Get a probability estimate for a value
   *
   * @param data the value to estimate the probability of
   * @param given the new value that data is conditional upon 
   * @return the estimated probability of the supplied value
   */
  public double getProbability(double data, double given) {

    return getEstimator(given).getProbability(data);
  }

  /**
   * Display a representation of this estimator
   */
  public String toString() {

    String result = "KK Conditional Estimator. " 
      + m_NumValues + " Normal Kernels:\n"
      + "StandardDev = " + Utils.doubleToString(m_StandardDev,4,2) 
      + "  \nMeans =";
    for(int i = 0; i < m_NumValues; i++) {
      result += " (" + m_Values[i] + ", " + m_CondValues[i] + ")";
      if (!m_AllWeightsOne) {
	  result += "w=" + m_Weights[i];
      }
    }
    return result;
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 1.8 $");
  }

  /**
   * Main method for testing this class. Creates some random points
   * in the range 0 - 100, 
   * and prints out a distribution conditional on some value
   *
   * @param argv should contain: seed conditional_value numpoints
   */
  public static void main(String [] argv) {

    try {
      int seed = 42;
      if (argv.length > 0) {
	seed = Integer.parseInt(argv[0]);
      }
      KKConditionalEstimator newEst = new KKConditionalEstimator(0.1);

      // Create 100 random points and add them
      Random r = new Random(seed);
      
      int numPoints = 50;
      if (argv.length > 2) {
	numPoints = Integer.parseInt(argv[2]);
      }
      for(int i = 0; i < numPoints; i++) {
	int x = Math.abs(r.nextInt()%100);
	int y = Math.abs(r.nextInt()%100);
	System.out.println("# " + x + "  " + y);
	newEst.addValue(x, y, 1);
      }
      //    System.out.println(newEst);
      int cond;
      if (argv.length > 1) {
	cond = Integer.parseInt(argv[1]);
      } else {
	cond = Math.abs(r.nextInt()%100);
      }
      System.out.println("## Conditional = " + cond);
      Estimator result = newEst.getEstimator(cond);
      for(int i = 0; i <= 100; i+= 5) {
	System.out.println(" " + i + "  " + result.getProbability(i));
      }
    } catch (Exception e) {
      System.out.println(e.getMessage());
    }
  }
}
