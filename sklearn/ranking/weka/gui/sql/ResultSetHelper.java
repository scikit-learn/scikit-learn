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
 * ResultSetHelper.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */


package weka.gui.sql;

import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.Types;
import java.util.Vector;

/**
 * Represents an extended JTable, containing a table model based on a ResultSet
 * and the corresponding query.
 *
 * @author      FracPete (fracpete at waikato dot ac dot nz)
 * @version     $Revision: 6429 $
 */
public class ResultSetHelper {
  
  /** the resultset to work on. */
  protected ResultSet m_ResultSet;

  /** whether we initialized. */
  protected boolean m_Initialized = false;

  /** the maximum number of rows to retrieve. */
  protected int m_MaxRows = 0;

  /** the number of columns. */
  protected int m_ColumnCount = 0;

  /** the number of rows. */
  protected int m_RowCount = 0;

  /** the column names. */
  protected String[] m_ColumnNames = null;

  /** whether a column is numeric. */
  protected boolean[] m_NumericColumns = null;

  /** the class for each column. */
  protected Class[] m_ColumnClasses = null;

  /**
   * initializes the helper, with unlimited number of rows.
   * 
   * @param rs        the resultset to work on
   */
  public ResultSetHelper(ResultSet rs) {
    this(rs, 0);
  }

  /**
   * initializes the helper, with the given maximum number of rows (less than
   * 1 means unlimited).
   * 
   * @param rs        the resultset to work on
   * @param max       the maximum number of rows to retrieve
   */
  public ResultSetHelper(ResultSet rs, int max) {
    super();

    m_ResultSet = rs;
    m_MaxRows   = max;
  }

  /**
   * initializes, i.e. reads the data, etc.
   */
  protected void initialize() {
    ResultSetMetaData     meta;
    int                   i;
    
    if (m_Initialized)
      return;
    
    try {
      meta = m_ResultSet.getMetaData();

      // columns names
      m_ColumnNames = new String[meta.getColumnCount()];
      for (i = 1; i <= meta.getColumnCount(); i++)
        m_ColumnNames[i - 1] = meta.getColumnName(i);
      
      // numeric columns
      m_NumericColumns = new boolean[meta.getColumnCount()];
      for (i = 1; i <= meta.getColumnCount(); i++)
        m_NumericColumns[i - 1] = typeIsNumeric(meta.getColumnType(i));

      // column classes
      m_ColumnClasses = new Class[meta.getColumnCount()];
      for (i = 1; i <= meta.getColumnCount(); i++) {
	try {
	  m_ColumnClasses[i - 1] = typeToClass(meta.getColumnType(i));
	}
	catch (Exception ex) {
	  m_ColumnClasses[i - 1] = String.class;
	}
      }

      // dimensions
      m_ColumnCount = meta.getColumnCount();

      // if the JDBC driver doesn't support scrolling we can't determine
      // the row count here
      if (m_ResultSet.getType() == ResultSet.TYPE_FORWARD_ONLY) {
	m_RowCount = -1;
      }
      else {
	m_RowCount = 0;
	m_ResultSet.first();
	if (m_MaxRows > 0) {
	  try {
	    m_ResultSet.absolute(m_MaxRows);
	    m_RowCount = m_ResultSet.getRow();
	  }
	  catch (Exception ex) {
	    // ignore it
	  }
	}
	else {
	  m_ResultSet.last();
	  m_RowCount = m_ResultSet.getRow();
	}

	// sometimes, e.g. with a "desc <table>", we can't use absolute(int)
	// and getRow()???
	try {
	  if ( (m_RowCount == 0) && (m_ResultSet.first()) ) {
	    m_RowCount = 1;
	    while (m_ResultSet.next()) {
	      m_RowCount++;
	      if (m_ResultSet.getRow() == m_MaxRows)
		break;
	    };
	  }
	}
	catch (Exception e) {
	  // ignore it
	}
      }

      m_Initialized = true;
    }
    catch (Exception ex) {
      // ignore it
    }
  }

  /**
   * the underlying resultset.
   * 
   * @return		the resultset
   */
  public ResultSet getResultSet() {
    return m_ResultSet;
  }

  /**
   * returns the number of columns in the resultset.
   * 
   * @return		the number of columns
   */
  public int getColumnCount() {
    initialize();

    return m_ColumnCount;
  }

  /**
   * returns the number of rows in the resultset. If -1 then the number of
   * rows couldn't be determined, i.e., the cursors aren't scrollable.
   * 
   * @return		the number of rows, -1 if it wasn't possible to 
   * 			determine
   */
  public int getRowCount() {
    initialize();

    return m_RowCount;
  }

  /**
   * returns an array with the names of the columns in the resultset.
   * 
   * @return		the column names
   */
  public String[] getColumnNames() {
    initialize();

    return m_ColumnNames;
  }

  /**
   * returns an array that indicates whether a column is numeric or nor.
   * 
   * @return		the numeric columns
   */
  public boolean[] getNumericColumns() {
    initialize();

    return m_NumericColumns;
  }

  /**
   * returns the classes for the columns.
   * 
   * @return		the column classes
   */
  public Class[] getColumnClasses() {
    initialize();

    return m_ColumnClasses;
  }

  /**
   * whether a limit on the rows to retrieve was set.
   * 
   * @return		true if there's a limit
   */
  public boolean hasMaxRows() {
    return (m_MaxRows > 0);
  }

  /**
   * the maximum number of rows to retrieve, less than 1 means unlimited.
   * 
   * @return		the maximum number of rows
   */
  public int getMaxRows() {
    return m_MaxRows;
  }

  /**
   * returns an 2-dimensional array with the content of the resultset, the first
   * dimension is the row, the second the column (i.e., getCells()[y][x]).
   * Note: the data is not cached! It is always retrieved anew.
   * 
   * @return		the data
   */
  public Object[][] getCells() {
    int			i;
    int			n;
    Vector<Object[]>	result;
    Object[]		row;
    int			rowCount;
    boolean		proceed;
    
    initialize();

    result = new Vector<Object[]>();

    try {
      
      // do know the number of rows?
      rowCount = getRowCount();
      if (rowCount == -1) {
	rowCount = getMaxRows();
	proceed  = m_ResultSet.next();
      }
      else {
	proceed = m_ResultSet.first();
      }
      
      if (proceed) {
	i = 0;
	while (true) {
	  row = new Object[getColumnCount()];
	  result.add(row);

	  for (n = 0; n < getColumnCount(); n++) {
	    try {
	      // to get around byte arrays when using getObject(int)
	      if (getColumnClasses()[n] == String.class)
		row[n] = m_ResultSet.getString(n + 1);
	      else
		row[n] = m_ResultSet.getObject(n + 1);
	    }
	    catch (Exception e) {
	      row[n] = null;
	    }
	  }

	  // get next row, if possible
	  if (i == rowCount - 1) {
	    break;
	  }
	  else {
	    // no more rows -> exit
	    if (!m_ResultSet.next())
	      break;
	  }

	  i++;
	}
      }
    }
    catch (Exception e) {
      e.printStackTrace();
    }

    return result.toArray(new Object[result.size()][getColumnCount()]);
  }

  /**
   * Returns the class associated with a SQL type.
   *
   * @param type the SQL type
   * @return the Java class corresponding with the type
   */
  public static Class typeToClass(int type) {
    Class     result;
    
    switch (type) {
      case Types.BIGINT :
        result = Long.class;
        break;
      case Types.BINARY:
        result = String.class;
        break;
      case Types.BIT:
        result = Boolean.class;
        break;
      case Types.CHAR:
        result = Character.class;
        break;
      case Types.DATE:
        result = java.sql.Date.class;
        break;
      case Types.DECIMAL:
        result = Double.class;
        break;
      case Types.DOUBLE:
        result = Double.class;
        break;
      case Types.FLOAT:
        result = Float.class;
        break;
      case Types.INTEGER:
        result = Integer.class;
        break;
      case Types.LONGVARBINARY:
        result = String.class;
        break;
      case Types.LONGVARCHAR:
        result = String.class;
        break;
      case Types.NULL:
        result = String.class;
        break;
      case Types.NUMERIC:
        result = Double.class;
        break;
      case Types.OTHER:
        result = String.class;
        break;
      case Types.REAL:
        result = Double.class;
        break;
      case Types.SMALLINT:
        result = Short.class;
        break;
      case Types.TIME:
        result = java.sql.Time.class;
        break;
      case Types.TIMESTAMP:
        result = java.sql.Timestamp.class;
        break;
      case Types.TINYINT:
        result = Short.class;
        break;
      case Types.VARBINARY:
        result = String.class;
        break;
      case Types.VARCHAR:
        result = String.class;
        break;
      default:
        result = null;
    }

    return result;
  }

  /**
   * returns whether the SQL type is numeric (and therefore the justification
   * should be right).
   * 
   * @param type      the SQL type
   * @return          whether the given type is numeric
   */
  public static boolean typeIsNumeric(int type) {
    boolean     result;
    
    switch (type) {
      case Types.BIGINT :
        result = true;
        break;
      case Types.BINARY:
        result = false;
      case Types.BIT:
        result = false;
        break;
      case Types.CHAR:
        result = false;
        break;
      case Types.DATE:
        result = false;
        break;
      case Types.DECIMAL:
        result = true;
        break;
      case Types.DOUBLE:
        result = true;
        break;
      case Types.FLOAT:
        result = true;
        break;
      case Types.INTEGER:
        result = true;
        break;
      case Types.LONGVARBINARY:
        result = false;
        break;
      case Types.LONGVARCHAR:
        result = false;
        break;
      case Types.NULL:
        result = false;
        break;
      case Types.NUMERIC:
        result = true;
        break;
      case Types.OTHER:
        result = false;
        break;
      case Types.REAL:
        result = true;
        break;
      case Types.SMALLINT:
        result = true;
        break;
      case Types.TIME:
        result = false;
        break;
      case Types.TIMESTAMP:
        result = true;
        break;
      case Types.TINYINT:
        result = true;
        break;
      case Types.VARBINARY:
        result = false;
        break;
      case Types.VARCHAR:
        result = false;
        break;
      default:
        result = false;
    }

    return result;
  }
}
