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
 *    StartSetHandler.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.attributeSelection;

/** 
 * Interface for search methods capable of doing something sensible
 * given a starting set of attributes.
 *
 * @author Mark Hall (mhall@cs.waikato.ac.nz)
 * @version $Revision: 1.8 $
 */
public interface StartSetHandler {

  /**
   * Sets a starting set of attributes for the search. It is the
   * search method's responsibility to report this start set (if any)
   * in its toString() method.
   * @param startSet a string containing a list of attributes (and or ranges),
   * eg. 1,2,6,10-15.
   * @exception Exception if start set can't be set.
   */
  void setStartSet (String startSet) throws Exception;

  /**
   * Returns a list of attributes (and or attribute ranges) as a String
   * @return a list of attributes (and or attribute ranges)
   */
  String getStartSet ();
}
