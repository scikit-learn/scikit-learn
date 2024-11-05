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
 * MultiInstanceClassifier.java
 * Copyright (C) 2006 University of Waikato, Hamilton, New Zealand
 */

package weka.core;

/**
 * Multi-Instance classifiers can specify an additional Capabilities object
 * for the data in the relational attribute, since the format of multi-instance
 * data is fixed to "bag/NOMINAL,data/RELATIONAL,class".
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 5953 $
 */
public interface MultiInstanceCapabilitiesHandler
  extends CapabilitiesHandler {

  /**
   * Returns the capabilities of this multi-instance classifier for the
   * relational data (i.e., the bags).
   *
   * @return            the capabilities of this object
   * @see               Capabilities
   */
  public Capabilities getMultiInstanceCapabilities();
}
