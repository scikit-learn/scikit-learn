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
 *    UnivariateDensityEstimator.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;

/**
 * Interface that can be implemented by simple weighted univariate
 * density estimators.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 5680 $
 */
public interface UnivariateDensityEstimator {

  /**
   * Adds a value to the density estimator.
   *
   * @param value the value to add
   * @param weight the weight of the value
   */
  void addValue(double value, double weight);

  /**
   * Returns the natural logarithm of the density estimate at the given
   * point.
   *
   * @param value the value at which to evaluate
   * @return the natural logarithm of the density estimate at the given
   * value
   */
  double logDensity(double value);
}