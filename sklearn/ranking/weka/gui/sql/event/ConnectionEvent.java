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
 * ConnectionEvent.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.sql.event;

import weka.gui.sql.DbUtils;

import java.util.EventObject;

/**
 * An event that is generated when a connection is established or dropped.
 *
 * @see         ConnectionListener
 * @author      FracPete (fracpete at waikato dot ac dot nz)
 * @version     $Revision: 1.2 $
 */
public class ConnectionEvent
  extends EventObject {

  /** for serialization */
  private static final long serialVersionUID = 5420308930427835037L;
  
  /** it was a connect try */
  public final static int CONNECT = 0;

  /** it was a disconnect */
  public final static int DISCONNECT = 1;
  
  /** the type of event, CONNECT or DISCONNECT */
  protected int m_Type;
  
  /** the databaseutils instance reponsible for the connection */
  protected DbUtils m_DbUtils;

  /** a possible exception that occurred if not successful */
  protected Exception m_Exception;
  
  /**
   * constructs the event
   * @param source        the source that generated this event
   * @param type          whether CONNECT or DISCONNECT happened
   * @param utils         the DatabaseUtils isntance responsible for the
   *                      connection
   */
  public ConnectionEvent(Object source, int type, DbUtils utils) {
    this(source, type, utils, null);
  }
  
  /**
   * constructs the event
   * @param source        the source that generated this event
   * @param type          whether CONNECT or DISCONNECT happened
   * @param utils         the DatabaseUtils isntance responsible for the
   *                      connection
   * @param ex            a possible exception, if not successful
   */
  public ConnectionEvent(Object source, int type, DbUtils utils, Exception ex) {
    super(source);
    
    m_Type      = type;
    m_DbUtils   = utils;
    m_Exception = ex;
  }
  
  /**
   * returns the type of this event, CONNECT or DISCONNECT
   * @return          the type of this event
   * @see             #CONNECT
   * @see             #DISCONNECT
   */
  public int getType() {
    return m_Type;
  }
  
  /**
   * whether an exception happened and is stored
   * @return          whether an exception happened
   */
  public boolean failed() {
    return (getException() != null);
  }
  
  /**
   * returns whether the connection is still open.
   * @return        whether the connection is still open
   */
  public boolean isConnected() {
    return m_DbUtils.isConnected();
  }

  /**
   * returns the stored exception, if any (can be NULL)
   */
  public Exception getException() {
    return m_Exception;
  }
  
  /**
   * returns the DbUtils instance that is responsible for the
   * connect/disconnect.
   * @return        the responsible DbUtils instance
   */
  public DbUtils getDbUtils() {
    return m_DbUtils;
  }

  /**
   * returns the event in a string representation
   * @return        the event in a string representation
   */
  public String toString() {
    String        result;

    result  = super.toString();
    result  = result.substring(0, result.length() - 1);  // remove "]"
    result +=   ",url=" + m_DbUtils.getDatabaseURL() 
              + ",user=" + m_DbUtils.getUsername()
              + ",password=" + m_DbUtils.getPassword().replaceAll(".", "*")
              + ",connected=" + isConnected() 
              + ",exception=" + getException()
              + "]";

    return result;
  }
}
