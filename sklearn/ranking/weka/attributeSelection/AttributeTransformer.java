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
 *    AttributeTransformer.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.attributeSelection;

import weka.core.Instance;
import weka.core.Instances;

/** 
 * Abstract attribute transformer. Transforms the dataset.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 4613 $
 */
public interface AttributeTransformer {
    // ===============
    // Public methods.
    // ===============

  /**
   * Returns just the header for the transformed data (ie. an empty
   * set of instances. This is so that AttributeSelection can
   * determine the structure of the transformed data without actually
   * having to get all the transformed data through getTransformedData().
   * @return the header of the transformed data.
   * @exception Exception if the header of the transformed data can't
   * be determined.
   */
  Instances transformedHeader() throws Exception;

  /**
   * Transform the supplied data set (assumed to be the same format
   * as the training data)
   * @return A set of instances representing the transformed data
   * @exception Exception if the attribute could not be evaluated
   */
  Instances transformedData(Instances data) throws Exception;

  /**
   * Transforms an instance in the format of the original data to the
   * transformed space
   * @return a transformed instance
   * @exception Exception if the instance could not be transformed
   */
  Instance convertInstance(Instance instance) throws Exception;
}
