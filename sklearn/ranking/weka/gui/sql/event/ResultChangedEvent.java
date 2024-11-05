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
 * ResultChangedEvent.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.sql.event;

import java.util.EventObject;

/**
 * An event that is generated when a different Result is activated in the
 * ResultPanel.
 *
 * @see         ResultChangedListener
 * @author      FracPete (fracpete at waikato dot ac dot nz)
 * @version     $Revision: 1.2 $
 */
public class ResultChangedEvent
  extends EventObject {

  /** for serialization */
  private static final long serialVersionUID = 36042516077236111L;
  
  /** the query that is associated with the active result table */
  protected String m_Query;

  /** the connect string with which the query was run */
  protected String m_URL;

  /** the user that was used to connect to the DB */
  protected String m_User;

  /** the password that was used to connect to the DB */
  protected String m_Password;

  /**
   * constructs the event
   * @param source        the source that generated this event
   * @param url           the current database url
   * @param user          the current user
   * @param pw            the current password
   * @param query         the current query
   */
  public ResultChangedEvent(Object source, 
                            String url,
                            String user,
                            String pw,
                            String query ) {
    super(source);

    m_URL      = url;
    m_User     = user;
    m_Password = pw;
    m_Query    = query;
  }

  /**
   * returns the database URL that produced the table model
   */
  public String getURL() {
    return m_URL;
  }

  /**
   * returns the user that produced the table model
   */
  public String getUser() {
    return m_User;
  }

  /**
   * returns the password that produced the table model
   */
  public String getPassword() {
    return m_Password;
  }

  /**
   * returns the query that was executed
   */
  public String getQuery() {
    return m_Query;
  }

  /**
   * returns the event in a string representation
   * @return        the event in a string representation
   */
  public String toString() {
    String        result;

    result  = super.toString();
    result  = result.substring(0, result.length() - 1);  // remove "]"
    result +=   ",url=" + getURL() 
              + ",user=" + getUser()
              + ",password=" + getPassword().replaceAll(".", "*")
              + ",query=" + getQuery()
              + "]";

    return result;
  }
}
