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
 *    NNConditionalEstimator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;

import java.util.Random;
import java.util.Vector;

import weka.core.matrix.Matrix;
import weka.core.RevisionUtils;
import weka.core.Utils;

/** 
 * Conditional probability estimator for a numeric domain conditional upon
 * a numeric domain (using Mahalanobis distance).
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.8 $
 */
public class NNConditionalEstimator implements ConditionalEstimator {

  /** Vector containing all of the values seen */
  private Vector m_Values = new Vector();

  /** Vector containing all of the conditioning values seen */
  private Vector m_CondValues = new Vector();

  /** Vector containing the associated weights */
  private Vector m_Weights = new Vector();

  /** The sum of the weights so far */
  private double m_SumOfWeights;

  /** Current Conditional mean */
  private double m_CondMean;

  /** Current Values mean */
  private double m_ValueMean;

  /** Current covariance matrix */
  private Matrix m_Covariance;

  /** Whether we can optimise the kernel summation */
  private boolean m_AllWeightsOne = true;

  /** 2 * PI */
  private static double TWO_PI = 2 * Math.PI;
  
  // ===============
  // Private methods
  // ===============

  /**
   * Execute a binary search to locate the nearest data value
   *
   * @param key the data value to locate
   * @param secondaryKey the data value to locate
   * @return the index of the nearest data value
   */
  private int findNearestPair(double key, double secondaryKey) {
    
    int low = 0; 
    int high = m_CondValues.size();
    int middle = 0;
    while (low < high) {
      middle = (low + high) / 2;
      double current = ((Double)m_CondValues.elementAt(middle)).doubleValue();
      if (current == key) {
	double secondary = ((Double)m_Values.elementAt(middle)).doubleValue();
	if (secondary == secondaryKey) {
	  return middle;
	}
	if (secondary > secondaryKey) {
	  high = middle;
	} else if (secondary < secondaryKey) {
	  low = middle + 1;
	}
      }
      if (current > key) {
	high = middle;
      } else if (current < key) {
	low = middle + 1;
      }
    }
    return low;
  }

  /** Calculate covariance and value means */
  private void calculateCovariance() {
    
    double sumValues = 0, sumConds = 0;
    for(int i = 0; i < m_Values.size(); i++) {
      sumValues += ((Double)m_Values.elementAt(i)).doubleValue()
	* ((Double)m_Weights.elementAt(i)).doubleValue();
      sumConds += ((Double)m_CondValues.elementAt(i)).doubleValue()
	* ((Double)m_Weights.elementAt(i)).doubleValue();
    }
    m_ValueMean = sumValues / m_SumOfWeights;
    m_CondMean = sumConds / m_SumOfWeights;
    double c00 = 0, c01 = 0, c10 = 0, c11 = 0;
    for(int i = 0; i < m_Values.size(); i++) {
      double x = ((Double)m_Values.elementAt(i)).doubleValue();
      double y = ((Double)m_CondValues.elementAt(i)).doubleValue();
      double weight = ((Double)m_Weights.elementAt(i)).doubleValue();
      c00 += (x - m_ValueMean) * (x - m_ValueMean) * weight;
      c01 += (x - m_ValueMean) * (y - m_CondMean) * weight;
      c11 += (y - m_CondMean) * (y - m_CondMean) * weight;
    }
    c00 /= (m_SumOfWeights - 1.0);
    c01 /= (m_SumOfWeights - 1.0);
    c10 = c01;
    c11 /= (m_SumOfWeights - 1.0);
    m_Covariance = new Matrix(2, 2);
    m_Covariance.set(0, 0, c00);
    m_Covariance.set(0, 1, c01);
    m_Covariance.set(1, 0, c10);
    m_Covariance.set(1, 1, c11);
  }

  /**
   * Returns value for normal kernel
   *
   * @param x the argument to the kernel function
   * @param variance the variance
   * @return the value for a normal kernel
   */
  private double normalKernel(double x, double variance) {
    
    return Math.exp(-x * x / (2 * variance)) / Math.sqrt(variance * TWO_PI);
  }
  
  /**
   * Add a new data value to the current estimator.
   *
   * @param data the new data value 
   * @param given the new value that data is conditional upon 
   * @param weight the weight assigned to the data value 
   */
  public void addValue(double data, double given, double weight) {
    
    int insertIndex = findNearestPair(given, data);
    if ((m_Values.size() <= insertIndex)
	|| (((Double)m_CondValues.elementAt(insertIndex)).doubleValue()
	    != given)
	|| (((Double)m_Values.elementAt(insertIndex)).doubleValue()
	    != data)) {
      m_CondValues.insertElementAt(new Double(given), insertIndex);
      m_Values.insertElementAt(new Double(data), insertIndex);
      m_Weights.insertElementAt(new Double(weight), insertIndex);
      if (weight != 1) {
	m_AllWeightsOne = false;
      }
    } else {
      double newWeight = ((Double)m_Weights.elementAt(insertIndex))
	.doubleValue();
      newWeight += weight;
      m_Weights.setElementAt(new Double(newWeight), insertIndex);
      m_AllWeightsOne = false;      
    }
    m_SumOfWeights += weight;
    // Invalidate any previously calculated covariance matrix
    m_Covariance = null;
  }

  /**
   * Get a probability estimator for a value
   *
   * @param given the new value that data is conditional upon 
   * @return the estimator for the supplied value given the condition
   */
  public Estimator getEstimator(double given) {
    
    if (m_Covariance == null) {
      calculateCovariance();
    }
    Estimator result = new MahalanobisEstimator(m_Covariance,
						given - m_CondMean,
						m_ValueMean);
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

  /** Display a representation of this estimator */
  public String toString() {
    
    if (m_Covariance == null) {
      calculateCovariance();
    }
    String result = "NN Conditional Estimator. "
      + m_CondValues.size() 
      + " data points.  Mean = " + Utils.doubleToString(m_ValueMean, 4, 2)
      + "  Conditional mean = " + Utils.doubleToString(m_CondMean, 4, 2);
    result += "  Covariance Matrix: \n" + m_Covariance;
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
   * Main method for testing this class.
   *
   * @param argv should contain a sequence of numeric values
   */
  public static void main(String [] argv) {
    
    try {
      int seed = 42;
      if (argv.length > 0) {
	seed = Integer.parseInt(argv[0]);
      }
      NNConditionalEstimator newEst = new NNConditionalEstimator();

      // Create 100 random points and add them
      Random r = new Random(seed);
      
      int numPoints = 50;
      if (argv.length > 2) {
	numPoints = Integer.parseInt(argv[2]);
      }
      for(int i = 0; i < numPoints; i++) {
	int x = Math.abs(r.nextInt() % 100);
	int y = Math.abs(r.nextInt() % 100);
	System.out.println("# " + x + "  " + y);
	newEst.addValue(x, y, 1);
      }
      //    System.out.println(newEst);
      int cond;
      if (argv.length > 1) {
	cond = Integer.parseInt(argv[1]);
      }
      else cond = Math.abs(r.nextInt() % 100);
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
