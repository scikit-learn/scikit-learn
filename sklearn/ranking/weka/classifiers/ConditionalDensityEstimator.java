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
 *    ConditionalDensityEstimator.java
 *    Copyright (C) 2009 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers;

import weka.core.Instance;

/**
 * Interface for numeric prediction schemes that can output conditional
 * density estimates.
 *
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 6041 $
 */
public interface ConditionalDensityEstimator {

  /**
   * Returns natural logarithm of density estimate for given value based on given instance.
   *
   * @param instance the instance to make the prediction for.
   * @param value the value to make the prediction for.
   * @return the natural logarithm of the density estimate
   * @exception Exception if the density cannot be computed
   */
  public double logDensity(Instance instance, double value) throws Exception;
}

