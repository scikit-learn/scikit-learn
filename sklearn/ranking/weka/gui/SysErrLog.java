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
 *    SysErrLog.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui;

import java.text.SimpleDateFormat;
import java.util.Date;

/** 
 * This Logger just sends messages to System.err.
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @version $Revision: 1.4 $
 */
public class SysErrLog implements Logger {

  /**
   * Gets a string containing current date and time.
   *
   * @return a string containing the date and time.
   */
  protected static String getTimestamp() {

    return (new SimpleDateFormat("yyyy.MM.dd hh:mm:ss")).format(new Date());
  }

  /**
   * Sends the supplied message to the log area. The current timestamp will
   * be prepended.
   *
   * @param message a value of type 'String'
   */
  public void logMessage(String message) {
    
    System.err.println("LOG " + SysErrLog.getTimestamp() + ": "
		       + message);
  }

  /**
   * Sends the supplied message to the status line.
   *
   * @param message the status message
   */
  public void statusMessage(String message) {

    System.err.println("STATUS: " + message);
  }
}
