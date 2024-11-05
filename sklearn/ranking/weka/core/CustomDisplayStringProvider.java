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

/**
 * CustomDisplayStringProvider.java
 * Copyright (C) 2010 University of Waikato, Hamilton, New Zealand
 */
package weka.core;

/**
 * For classes that do not implement the OptionHandler interface and want to 
 * provide a custom display string in the GenericObjectEditor, which is more 
 * descriptive than the class name.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6287 $
 */
public interface CustomDisplayStringProvider {

  /**
   * Returns the custom display string.
   * 
   * @return		the string
   */
  public String toDisplay();
}
