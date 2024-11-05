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
 * Logger.java
 * Copyright (C) 2008 University of Waikato, Hamilton, New Zealand
 */

package weka.core.logging;

import weka.core.RevisionHandler;
import weka.core.Utils;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.SimpleDateFormat;
import java.util.Properties;

/**
 * Abstract superclass for all loggers.
 * 
 * @author  fracpete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 4716 $
 */
public abstract class Logger
  implements RevisionHandler {

  /** the properties file. */
  public final static String PROPERTIES_FILE = "weka/core/logging/Logging.props";
  
  /**
   * The logging level.
   * 
   * @author  fracpete (fracpete at waikato dot ac dot nz)
   * @version $Revision: 4716 $
   */
  public enum Level {
    /** logs all messages. */
    ALL(0),
    /** FINEST level. */
    FINEST(1),
    /** FINEST level. */
    FINER(2),
    /** FINER level. */
    FINE(3),
    /** FINE level. */
    INFO(4),
    /** WARNING level. */
    WARNING(5),
    /** SEVERE level. */
    SEVERE(6),
    /** turns logging off. */
    OFF(10);
    
    /** the order of the level. */
    private int m_Order;
    
    /**
     * Initializes the level.
     * 
     * @param order	the order of the level
     */
    private Level(int order) {
      m_Order = order;
    }
    
    /**
     * Returns the order of this level.
     * 
     * @return		the order
     */
    public int getOrder() {
      return m_Order;
    }
  }
  
  /** the minimum level of log events to have in order to end up in the log. */
  protected Level m_MinLevel;

  /** the singleton instance of the logger. */
  protected static Logger m_Singleton;
  
  /** the properties file. */
  protected static Properties m_Properties;

  /** for formatting the dates. */
  protected static SimpleDateFormat m_DateFormat;
  
  static {
    try {
      m_Properties = Utils.readProperties(PROPERTIES_FILE);
    }
    catch (Exception e) {
      System.err.println(
	  "Error reading the logging properties '" + PROPERTIES_FILE + "': " + e);
      m_Properties = new Properties();
    }
  }
  
  /**
   * Initializes the logger.
   */
  public Logger() {
    super();
    
    initialize();
  }
  
  /**
   * Initializes the logger.
   */
  protected void initialize() {
    m_MinLevel = Level.valueOf(m_Properties.getProperty("MinLevel", "INFO"));
  }
  
  /**
   * Returns the minimum level log messages must have in order to appear in
   * the log.
   * 
   * @return		the level
   */
  public Level getMinLevel() {
    return m_MinLevel;
  }
  
  /**
   * Returns the location the logging happened.
   * 
   * @return		the classname (= [0]), the method (= [1]) and the
   * 			line number (= [2]) that generated the logging event
   */
  protected static String[] getLocation() {
    String[]		result;
    Throwable 		t;
    StackTraceElement[]	trace;
    int			i;
    
    result = new String[3];
    
    t = new Throwable();
    t.fillInStackTrace();
    trace = t.getStackTrace();

    for (i = 0; i < trace.length; i++) {
      // skip the Logger class
      if (trace[i].getClassName().equals(Logger.class.getName()))
	continue;
      
      if (trace[i].getClassName().equals(weka.gui.LogPanel.class.getName()))
        continue;
      
      // fill in result
      result[0] = trace[i].getClassName();
      result[1] = trace[i].getMethodName();
      result[2] = "" + trace[i].getLineNumber();
      break;
    }
    
    return result;
  }
  
  /**
   * Performs the actual logging. 
   * Actual logger implementations must override this method.
   * 
   * @param level	the level of the message
   * @param msg		the message to log
   * @param cls		the classname originating the log event
   * @param method	the method originating the log event
   * @param lineno	the line number originating the log event
   */
  protected abstract void doLog(Level level, String msg, String cls, String method, int lineno);
  
  /**
   * Returns the singleton instance of the logger.
   * 
   * @return		the logger instance
   */
  public static Logger getSingleton() {
    String	classname;
    
    if (m_Singleton == null) {
      // logger
      classname = m_Properties.getProperty("Logger", ConsoleLogger.class.getName());
      try {
	m_Singleton = (Logger) Class.forName(classname).newInstance();
      }
      catch (Exception e) {
	e.printStackTrace();
      }
      
      // date format
      m_DateFormat = new SimpleDateFormat(m_Properties.getProperty("DateFormat", "yyyy-MM-dd HH:mm:ss"));
    }
    
    return m_Singleton;
  }
  
  /**
   * Logs the given message under the given level.
   * 
   * @param level	the level of the message
   * @param msg		the message to log
   */
  public static void log(Level level, String msg) {
    Logger	logger;
    boolean	log;
    String[]	location;
    
    logger = getSingleton();
    if (logger == null)
      return;
    
    synchronized(logger) {
      log = false;
      if (logger.getMinLevel() == Level.ALL)
	log = true;
      else if (level.getOrder() >= logger.getMinLevel().getOrder())
	log = true;
      if (!log)
	return;
      location = getLocation();
      logger.doLog(level, msg, location[0], location[1], Integer.parseInt(location[2]));
    }
  }
  
  /**
   * Logs the given message under the given level.
   * 
   * @param level	the level of the message
   * @param t		the throwable to log
   */
  public static void log(Level level, Throwable t) {
    StringWriter	swriter;
    PrintWriter		pwriter;
    
    swriter = new StringWriter();
    pwriter = new PrintWriter(swriter);
    t.printStackTrace(pwriter);
    pwriter.close();
    
    log(level, swriter.toString());
  }
}
