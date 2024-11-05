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
 *    Associator.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.associations;

import weka.core.Instances;
import weka.core.Capabilities;

public interface Associator {

  /**
   * Generates an associator. Must initialize all fields of the associator
   * that are not being set via options (ie. multiple calls of buildAssociator
   * must always lead to the same result). Must not change the dataset
   * in any way.
   *
   * @param data set of instances serving as training data 
   * @exception Exception if the associator has not been 
   * generated successfully
   */
  void buildAssociations(Instances data) throws Exception;

  /** 
   * Returns the Capabilities of this associator. Derived associators have to
   * override this method to enable capabilities.
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  Capabilities getCapabilities();
}