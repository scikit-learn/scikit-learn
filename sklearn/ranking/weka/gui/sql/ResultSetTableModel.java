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
 * ResultSetTableModel.java
 * Copyright (C) 2005 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui.sql;

import java.sql.ResultSet;
import java.util.HashSet;

import javax.swing.event.TableModelListener;
import javax.swing.table.TableModel;

/**
* The model for an SQL ResultSet.
*
* @author     FracPete (fracpete at waikato dot ac dot nz)
* @version    $Revision: 1.3 $
*/
public class ResultSetTableModel implements TableModel {
  
  /** the listeners. */
  protected HashSet m_Listeners;
  
  /** the data. */
  protected Object[][] m_Data;
  
  /** for retrieving the data etc. */
  protected ResultSetHelper m_Helper;

  /**
   * initializes the model, retrieves all rows.
   * 
   * @param rs          the ResultSet to get the data from
   */
  public ResultSetTableModel(ResultSet rs) {
    this(rs, 0);
  }

  /**
   * initializes the model, retrieves only the given amount of rows (0 means
   * all).
   * 
   * @param rs          the ResultSet to get the data from
   * @param rows        the maximum number of rows to retrieve, 0 retrieves all
   */
  public ResultSetTableModel(ResultSet rs, int rows) {
    super();

    m_Listeners = new HashSet();
    m_Helper    = new ResultSetHelper(rs, rows);
    m_Data      = m_Helper.getCells();
  }

  /**
   * adds a listener to the list that is notified each time a change to data 
   * model occurs.
   * 
   * @param l		the listener to add
   */
  public void addTableModelListener(TableModelListener l) {
    m_Listeners.add(l);
  }

  /**
   * returns the most specific superclass for all the cell values in the 
   * column (always String).
   * 
   * @param columnIndex	the index of the column
   * @return		the class
   */
  public Class getColumnClass(int columnIndex) {
    Class       result;

    result = null;

    if (    (m_Helper.getColumnClasses() != null) 
         && (columnIndex >= 0) 
         && (columnIndex < getColumnCount()) ) {
      if (columnIndex == 0)
        result = Integer.class;
      else
        result = m_Helper.getColumnClasses()[columnIndex - 1];
   }

    return result;
  }

  /**
   * returns the number of columns in the model.
   * 
   * @return		the number of columns
   */
  public int getColumnCount() {
    return m_Helper.getColumnCount() + 1;
  }

  /**
   * returns the name of the column at columnIndex.
   * 
   * @param columnIndex	the index of the column
   * @return		the name
   */
  public String getColumnName(int columnIndex) {
    String         result;

    result = "";

    if (    (m_Helper.getColumnNames() != null) 
        && (columnIndex >= 0) 
        && (columnIndex < getColumnCount()) ) {
      if (columnIndex == 0)
        result = "Row";
      else
        result = m_Helper.getColumnNames()[columnIndex - 1];
    }

    return result;
  }

  /**
   * returns the number of rows in the model.
   * 
   * @return		the number of data rows
   */
  public int getRowCount() {
    return m_Data.length;
  }

  /**
   * returns the value for the cell at columnindex and rowIndex.
   * 
   * @param rowIndex	the row of the cell
   * @param columnIndex	the column of the cell
   * @return		the data value
   */
  public Object getValueAt(int rowIndex, int columnIndex) {
    Object            result;

    result = null;

    if (    (rowIndex >= 0) && (rowIndex < getRowCount())
         && (columnIndex >= 0) && (columnIndex < getColumnCount()) ) {
      if (columnIndex == 0)
        result = new Integer(rowIndex + 1);
      else
        result = m_Data[rowIndex][columnIndex - 1];
    }

    return result;
  }

  /**
   * checks whether the value of the cell is NULL.
   * 
   * @param rowIndex	the row of the cell
   * @param columnIndex	the column of the cell
   * @return		true if the cell value is NULL
   */
  public boolean isNullAt(int rowIndex, int columnIndex) {
    return (getValueAt(rowIndex, columnIndex) == null);
  }

  /**
   * returns whether the column at the given index is numeric.
   * 
   * @param columnIndex       the column to check
   * @return                  whether the column is numeric
   */
  public boolean isNumericAt(int columnIndex) {
    boolean         result;

    result = false;
    
    if ( (columnIndex >= 0) && (columnIndex < getColumnCount()) ) {
      if (columnIndex == 0) {
        result = true;
      }
      else {
        if (m_Helper.getNumericColumns() == null)
          result = false;
        else
          result = m_Helper.getNumericColumns()[columnIndex - 1];
      }
    }

    return result;
  }

  /**
   * returns true if the cell at rowindex and columnindexis editable.
   * 
   * @param rowIndex	the row of the cell
   * @param columnIndex	the column of the cell
   * @return		always false
   */
  public boolean isCellEditable(int rowIndex, int columnIndex) {
    return false;
  }

  /**
   * removes a listener from the list that is notified each time a change to
   * the data model occurs.
   * 
   * @param l		the listener to remove
   */
  public void removeTableModelListener(TableModelListener l) {
    m_Listeners.remove(l);
  }

  /**
   * sets the value in the cell at columnIndex and rowIndex to aValue.
   * Ignored.
   * 
   * @param aValue	the value to set - ignored
   * @param rowIndex	the row of the cell
   * @param columnIndex	the column of the cell
   */
  public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
    // ignore
  }

  /**
   * frees up the memory.
   * 
   * @throws Throwable	if something goes wrong
   */
  public void finalize() throws Throwable {
    try {
      m_Helper.getResultSet().close();
      m_Helper.getResultSet().getStatement().close();
      m_Helper = null;
    }
    catch (Exception e) {
      // ignored
    }

    m_Data = null;
    
    super.finalize();
  }
}
