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
 *    Prediction.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.classifiers.evaluation;

/**
 * Encapsulates a single evaluatable prediction: the predicted value plus the 
 * actual class value.
 *
 * @author Len Trigg (len@reeltwo.com)
 * @version $Revision: 5987 $
 */
public interface Prediction {

  /** 
   * Constant representing a missing value. This should have the same value
   * as weka.core.Instance.MISSING_VALUE 
   */
  double MISSING_VALUE 
    = weka.core.Utils.missingValue();

  /** 
   * Gets the weight assigned to this prediction. This is typically the weight
   * of the test instance the prediction was made for.
   *
   * @return the weight assigned to this prediction.
   */
  double weight();

  /** 
   * Gets the actual class value.
   *
   * @return the actual class value, or MISSING_VALUE if no
   * prediction was made.  
   */
  double actual();

  /**
   * Gets the predicted class value.
   *
   * @return the predicted class value, or MISSING_VALUE if no
   * prediction was made.  
   */
  double predicted();

}
