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
 *    UnassignedDatasetException.java
 *    Copyright (C) 2002 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core;

/**
 * Exception that is raised when trying to use something that has no
 * reference to a dataset, when one is required.
 *
 * @author Richard Kirkby (rkirkby@cs.waikato.ac.nz)
 * @version $Revision: 5953 $
 */
public class UnassignedDatasetException
  extends RuntimeException {

  /** for serialization */
  private static final long serialVersionUID = -9000116786626328854L;

  /**
   * Creates a new UnassignedDatasetException with no message.
   *
   */
  public UnassignedDatasetException() {

    super();
  }

  /**
   * Creates a new UnassignedDatasetException.
   *
   * @param message the reason for raising an exception.
   */
  public UnassignedDatasetException(String message) {

    super(message);
  }
}
