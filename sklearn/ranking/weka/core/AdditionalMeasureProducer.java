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
 *    AdditionalMeasureProducer.java
 *    Copyright (C) 2000 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

import java.util.*;

/** 
 * Interface to something that can produce measures other than those
 * calculated by evaluation modules. 
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 5953 $
 */
public interface AdditionalMeasureProducer {

  /**
   * Returns an enumeration of the measure names. Additional measures
   * must follow the naming convention of starting with "measure", eg.
   * double measureBlah()
   * @return an enumeration of the measure names
   */
  Enumeration enumerateMeasures();

  /**
   * Returns the value of the named measure
   * @param measureName the name of the measure to query for its value
   * @return the value of the named measure
   * @exception IllegalArgumentException if the named measure is not supported
   */
  double getMeasure(String measureName);
}
