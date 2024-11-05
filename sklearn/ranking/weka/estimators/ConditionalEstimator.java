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
 *    ConditionalEstimator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;

import weka.core.RevisionHandler;

 
/** 
 * Interface for conditional probability estimators. Example code: <p>
 *
 * <code> <pre>
 *   NNConditionalEstimator newEst = new NNConditionalEstimator();
 *
 *   // Create 50 random points and add them
 *   Random r = new Random(seed);
 *   for(int i = 0; i < 50; i++) {
 *     int x = Math.abs(r.nextInt() % 100);
 *     int y = Math.abs(r.nextInt() % 100);
 *     System.out.println("# " + x + "  " + y);
 *     newEst.addValue(x, y, 1);
 *   }
 *
 *   // Pick a random conditional value
 *   int cond = Math.abs(r.nextInt() % 100);
 *   System.out.println("## Conditional = " + cond);
 *
 *   // Print the probabilities conditional on that value
 *   Estimator result = newEst.getEstimator(cond);
 *   for(int i = 0; i <= 100; i+= 5) {
 *     System.out.println(" " + i + "  " + result.getProbability(i));
 *   }
 * </pre> </code>
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.7 $
 */
public interface ConditionalEstimator extends RevisionHandler {

  /**
   * Add a new data value to the current estimator.
   *
   * @param data the new data value 
   * @param given the new value that data is conditional upon 
   * @param weight the weight assigned to the data value 
   */
  void addValue(double data, double given, double weight);

  /**
   * Get a probability estimator for a value
   *
   * @param given the new value that data is conditional upon 
   * @return the estimator for the supplied value given the condition
   */
  Estimator getEstimator(double given);

  /**
   * Get a probability for a value conditional on another value
   *
   * @param data the value to estimate the probability of
   * @param given the new value that data is conditional upon 
   * @return the estimator for the supplied value given the condition
   */
  double getProbability(double data, double given);

}








