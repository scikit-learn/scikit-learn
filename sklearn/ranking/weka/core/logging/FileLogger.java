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
 * FileLogger.java
 * Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 */

package weka.core.logging;

import weka.core.RevisionUtils;
import weka.core.WekaPackageManager;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Date;

/**
 * A simple file logger, that just logs to a single file. Deletes the file
 * when an object gets instantiated.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6679 $
 */
public class FileLogger
  extends ConsoleLogger {

  /** the log file. */
  protected File m_LogFile;
  
  /** the line feed. */
  protected String m_LineFeed;
  
  /**
   * Initializes the logger.
   */
  protected void initialize() {
    super.initialize();

    // log file
    m_LogFile = getLogFile();
    // try to remove file
    try {
      if ((m_LogFile != null) && m_LogFile.exists())
	m_LogFile.delete();
    }
    catch (Exception e) {
      e.printStackTrace();
    }
    
    // the line feed
    m_LineFeed = System.getProperty("line.separator");
  }
  
  /**
   * Returns the log file to use.
   * 
   * @return		the log file
   */
  protected File getLogFile() {
    String	filename;
    File	result;
    
    filename = m_Properties.getProperty("LogFile", "%w/weka.log");
    filename = filename.replaceAll("%t", System.getProperty("java.io.tmpdir"));
    filename = filename.replaceAll("%h", System.getProperty("user.home"));
    filename = filename.replaceAll("%c", System.getProperty("user.dir"));
    filename = filename.replaceAll("%w", WekaPackageManager.WEKA_HOME.toString());
    filename = filename.replaceAll("%%", System.getProperty("%"));
    
    result = new File(filename);
    
    return result;
  }
  
  /**
   * Appends the given string to the log file (without new line!).
   * 
   * @param s		the string to append
   */
  protected void append(String s) {
    BufferedWriter	writer;
   
    if (m_LogFile == null)
      return;
    
    // append output to file
    try {
      writer = new BufferedWriter(new FileWriter(m_LogFile, true));
      writer.write(s);
      writer.flush();
      writer.close();
    }
    catch (Exception e) {
      // ignored
    }
  }

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
    // output to console
    super.doLog(level, msg, cls, method, lineno);
    
    // append output to file
    append(
	m_DateFormat.format(new Date()) + " " + cls + " " + method + m_LineFeed
	+ level + ": " + msg + m_LineFeed);
  }
  
  /**
   * Returns the revision string.
   * 
   * @return		the revision
   */
  public String getRevision() {
    return RevisionUtils.extract("$Revision: 6679 $");
  }
}
