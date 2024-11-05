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
 *    IncrementalEstimator.java
 *    Copyright (C) 2004 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.estimators;


/** 
 * Interface for an incremental probability estimators.<p>
 *
 * @author Gabi Schmidberger (gabi@cs.waikato.ac.nz)
 * @version $Revision: 1.2 $
 */
public interface IncrementalEstimator {

  /**
   * Add one value to the current estimator.
   *
   * @param data the new data value 
   * @param weight the weight assigned to the data value 
   */
  void addValue(double data, double weight);

}








