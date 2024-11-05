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
 *    UnsuppotedClassTypeException.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

/**
 * Exception that is raised by an object that is unable to process the
 * class type of the data it has been passed.
 *
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @version $Revision: 5953 $
 */
public class UnsupportedClassTypeException
  extends WekaException {

  /** for serialization */
  private static final long serialVersionUID = 5175741076972192151L;

  /**
   * Creates a new UnsupportedClassTypeException with no message.
   *
   */
  public UnsupportedClassTypeException() {

    super();
  }

  /**
   * Creates a new UnsupportedClassTypeException.
   *
   * @param message the reason for raising an exception.
   */
  public UnsupportedClassTypeException(String message) {

    super(message);
  }
}
