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
 *    Logger.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui;

/** 
 * Interface for objects that display log (permanent historical) and
 * status (transient) messages.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */
public interface Logger {

  /**
   * Sends the supplied message to the log area. These message will typically
   * have the current timestamp prepended, and be viewable as a history.
   *
   * @param message the log message
   */
  void logMessage(String message);
  
  /**
   * Sends the supplied message to the status line. These messages are
   * typically one-line status messages to inform the user of progress
   * during processing (i.e. it doesn't matter if the user doesn't happen
   * to look at each message)
   *
   * @param message the status message.
   */
  void statusMessage(String message);
  
}
