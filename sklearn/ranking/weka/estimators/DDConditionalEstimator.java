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
 *    DDConditionalEstimator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;

import weka.core.RevisionUtils;

 
/** 
 * Conditional probability estimator for a discrete domain conditional upon
 * a discrete domain.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.8 $
 */
public class DDConditionalEstimator implements ConditionalEstimator {

  /** Hold the sub-estimators */
  private DiscreteEstimator [] m_Estimators;

  /**
   * Constructor
   *
   * @param numSymbols the number of possible symbols (remember to include 0)
   * @param numCondSymbols the number of conditioning symbols 
   * @param laplace if true, sub-estimators will use laplace
   */
  public DDConditionalEstimator(int numSymbols, int numCondSymbols,
				boolean laplace) {
    
    m_Estimators = new DiscreteEstimator [numCondSymbols];
    for(int i = 0; i < numCondSymbols; i++) {
      m_Estimators[i] = new DiscreteEstimator(numSymbols, laplace);
    }
  }

  /**
   * Add a new data value to the current estimator.
   *
   * @param data the new data value 
   * @param given the new value that data is conditional upon 
   * @param weight the weight assigned to the data value 
   */
  public void addValue(double data, double given, double weight) {
    
    m_Estimators[(int)given].addValue(data, weight);
  }

  /**
   * Get a probability estimator for a value
   *
   * @param given the new value that data is conditional upon 
   * @return the estimator for the supplied value given the condition
   */
  public Estimator getEstimator(double given) {
    
    return m_Estimators[(int)given];
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
    
    String result = "DD Conditional Estimator. " 
      + m_Estimators.length + " sub-estimators:\n";
    for(int i = 0; i < m_Estimators.length; i++) {
      result += "Sub-estimator " + i + ": " + m_Estimators[i];
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
   * Main method for testing this class.
   *
   * @param argv should contain a sequence of pairs of integers which
   * will be treated as symbolic.
   */
  public static void main(String [] argv) {
    
    try {
      if (argv.length == 0) {
	System.out.println("Please specify a set of instances.");
	return;
      }
      int currentA = Integer.parseInt(argv[0]);
      int maxA = currentA;
      int currentB = Integer.parseInt(argv[1]);
      int maxB = currentB;
      for(int i = 2; i < argv.length - 1; i += 2) {
	currentA = Integer.parseInt(argv[i]);
	currentB = Integer.parseInt(argv[i + 1]);
	if (currentA > maxA) {
	  maxA = currentA;
	}
	if (currentB > maxB) {
	  maxB = currentB;
	}
      }
      DDConditionalEstimator newEst = new DDConditionalEstimator(maxA + 1,
								 maxB + 1,
								 true);
      for(int i = 0; i < argv.length - 1; i += 2) {
	currentA = Integer.parseInt(argv[i]);
	currentB = Integer.parseInt(argv[i + 1]);
	System.out.println(newEst);
	System.out.println("Prediction for " + currentA + '|' + currentB 
			   + " = "
			   + newEst.getProbability(currentA, currentB));
	newEst.addValue(currentA, currentB, 1);
      }
    } catch (Exception e) {
      System.out.println(e.getMessage());
    }
  }
}
