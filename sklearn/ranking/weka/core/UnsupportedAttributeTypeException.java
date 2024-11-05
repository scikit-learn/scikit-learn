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
 *    UnsuppotedAttributeTypeException.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

/**
 * Exception that is raised by an object that is unable to process some of the
 * attribute types it has been passed.
 *
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @version $Revision: 5953 $
 */
public class UnsupportedAttributeTypeException
  extends WekaException {

  /** for serialization */
  private static final long serialVersionUID = 2658987325328414838L;

  /**
   * Creates a new UnsupportedAttributeTypeException with no message.
   *
   */
  public UnsupportedAttributeTypeException() {

    super();
  }

  /**
   * Creates a new UnsupportedAttributeTypeException.
   *
   * @param message the reason for raising an exception.
   */
  public UnsupportedAttributeTypeException(String message) {

    super(message);
  }
}
