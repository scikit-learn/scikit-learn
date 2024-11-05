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
 *    IntervalEstimator.java
 *    Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers;

import weka.core.Instance;

/**
 * Interface for numeric prediction schemes that can output prediction
 * intervals.
 *
 * @author Kurt Driessens (kurtd@cs.waikato.ac.nz)
 * @version $Revision: 6041 $
 */
public interface IntervalEstimator {

  /**
   * Returns an N * 2 array, where N is the number of prediction
   * intervals. In each row, the first element contains the lower
   * boundary of the corresponding prediction interval and the second
   * element the upper boundary.
   *
   * @param inst the instance to make the prediction for.
   * @param confidenceLevel the percentage of cases that the interval should cover.
   * @return an array of prediction intervals
   * @exception Exception if the intervals can't be computed
   */
  double[][] predictIntervals(Instance inst, double confidenceLevel) throws Exception;
}

