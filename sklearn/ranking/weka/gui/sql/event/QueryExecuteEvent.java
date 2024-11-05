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
 * QueryExecuteEvent.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.sql.event;

import weka.gui.sql.DbUtils;

import java.sql.ResultSet;
import java.util.EventObject;

/**
 * An event that is generated when a query is executed.
 *
 * @see         QueryExecuteListener
 * @author      FracPete (fracpete at waikato dot ac dot nz)
 * @version     $Revision: 1.2 $
 */
public class QueryExecuteEvent
  extends EventObject {

  /** for serialization */
  private static final long serialVersionUID = -5556385019954730740L;
  
  /** the Db utils instance for the current DB connection  */
  protected DbUtils m_DbUtils;
  
  /** the query that was executed */
  protected String m_Query;

  /** the produced ResultSet, if any */
  protected ResultSet m_ResultSet;

  /** a possible exception, if the query failed */
  protected Exception m_Exception;

  /** the maximum number of rows to retrieve */
  protected int m_MaxRows;
  
  /**
   * constructs the event
   * @param source        the source that generated this event
   * @param utils         the DbUtils instance that connected to the DB
   * @param query         the query that is the basis for the resultset
   * @param rows          the maximum number of rows to retrieve (0 for all)
   * @param rs            the ResultSet that was produced (depending on the
   *                      type of SQL query it can also be NULL)
   * @param ex            in case an exception occurred
   */
  public QueryExecuteEvent( Object source, 
                            DbUtils utils,
                            String query, 
                            int rows,
                            ResultSet rs, 
                            Exception ex ) {
    super(source);

    m_DbUtils   = utils;
    m_Query     = query;
    m_MaxRows   = rows;
    m_ResultSet = rs;
    m_Exception = ex;
  }

  /**
   * returns the DbUtils instance that was executed the query
   */
  public DbUtils getDbUtils() {
    return m_DbUtils;
  }

  /**
   * returns the query that was executed
   */
  public String getQuery() {
    return m_Query;
  }

  /**
   * returns the maximum number of rows to retrieve. 0 means all.
   */
  public int getMaxRows() {
    return m_MaxRows;
  }

  /**
   * is TRUE in case the exception is not NULL, i.e. the query failed
   */
  public boolean failed() {
    return (m_Exception != null);
  }

  /**
   * whether a ResultSet was produced, e.g. DDL commands like delete, drop
   * or update do not produce one.
   */
  public boolean hasResult() {
    return (m_ResultSet != null);
  }

  /**
   * returns the resultset that was produced, can be null in case the query
   * failed
   */
  public ResultSet getResultSet() {
    return m_ResultSet;
  }

  /**
   * returns the exception, if one happened, otherwise NULL
   */
  public Exception getException() {
    return m_Exception;
  }

  /**
   * returns the event in a string representation
   * @return        the event in a string representation
   */
  public String toString() {
    String        result;

    result  = super.toString();
    result  = result.substring(0, result.length() - 1);  // remove "]"
    result +=   ",query=" + getQuery() 
              + ",maxrows=" + getMaxRows()
              + ",failed=" + failed()
              + ",exception=" + getException() + "]";

    return result;
  }
}
