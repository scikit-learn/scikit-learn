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
 *    WekaWrapper.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.beans;

/**
 * Interface to something that can wrap around a class of Weka
 * algorithms (classifiers, filters etc). Typically implemented
 * by a bean for handling classes of Weka algorithms. 
 *
 * @author <a href="mailto:mhall@cs.waikato.ac.nz">Mark Hall</a>
 * @version $Revision: 1.4 $
 * @since 1.0
 */
public interface WekaWrapper {

  /**
   * Set the algorithm.
   *
   * @param algorithm an <code>Object</code> value
   * @exception IllegalArgumentException if the supplied object is
   * not of the class of algorithms handled by this wrapper.
   */
  void setWrappedAlgorithm(Object algorithm);

  /**
   * Get the algorithm
   *
   * @return an <code>Object</code> value
   */
  Object getWrappedAlgorithm();
}
