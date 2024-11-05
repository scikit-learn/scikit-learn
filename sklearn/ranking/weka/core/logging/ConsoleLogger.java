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
 * ConsoleLogger.java
 * Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 */

package weka.core.logging;

import weka.core.RevisionUtils;

import java.util.Date;

/**
 * A simple logger that outputs the logging information in the console.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 4716 $
 */
public class ConsoleLogger
  extends Logger {

  /**
   * Performs the actual logging. 
   * 
   * @param level	the level of the message
   * @param msg		the message to log
   * @param cls		the classname originating the log event
   * @param method	the method originating the log event
   * @param lineno	the line number originating the log event
   */
  protected void doLog(Level level, String msg, String cls, String method, int lineno) {
    System.err.println(
	m_DateFormat.format(new Date()) + " " + cls + " " + method + "\n" 
	+ level + ": " + msg);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 4716 $");
  }
}
