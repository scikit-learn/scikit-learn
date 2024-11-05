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
 * SortedTableModel.java
 * Copyright (C) 2005-2010 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.gui;

import java.awt.event.InputEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collections;

import javax.swing.JTable;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumnModel;
import javax.swing.table.TableModel;

import weka.core.ClassDiscovery;

/**
 * Represents a TableModel with sorting functionality.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 * @version $Revision: 6275 $
 */

public class SortedTableModel
  extends AbstractTableModel
  implements TableModelListener {

  /** for serialization */
  static final long serialVersionUID = 4030907921461127548L;

  /**
   * Helper class for sorting the columns.
   */
  public static class SortContainer
    implements Comparable<SortContainer> {

    /** the value to sort. */
    protected Comparable m_Value;

    /** the index of the value. */
    protected int m_Index;

    /**
     * Initializes the container.
     *
     * @param value	the value to sort on
     * @param index	the original index
     */
    public SortContainer(Comparable value, int index) {
      super();

      m_Value = value;
      m_Index = index;
    }

    /**
     * Returns the value to sort on.
     *
     * @return		the value
     */
    public Comparable getValue() {
      return m_Value;
    }

    /**
     * Returns the original index of the item.
     *
     * @return		the index
     */
    public int getIndex() {
      return m_Index;
    }

    /**
     * Compares this object with the specified object for order.  Returns a
     * negative integer, zero, or a positive integer as this object is less
     * than, equal to, or greater than the specified object. Null is considered
     * smallest. If both values are null, then 0 is returned.
     *
     * @param o 	the object to be compared.
     * @return  	a negative integer, zero, or a positive integer as this object
     *			is less than, equal to, or greater than the specified object.
     * @throws ClassCastException 	if the specified object's type prevents it
     *         				from being compared to this object.
     */
    public int compareTo(SortContainer o) {
      if ((m_Value == null) || (o.getValue() == null)) {
	if (m_Value == o.getValue())
	  return 0;
	if (m_Value == null)
	  return -1;
	else
	  return +1;
      }
      else {
	return m_Value.compareTo(o.getValue());
      }
    }

    /**
     * Indicates whether some other object is "equal to" this one.
     *
     * @param obj	the reference object with which to compare.
     * @return		true if this object is the same as the obj argument;
     * 			false otherwise.
     * @throws ClassCastException 	if the specified object's type prevents it
     *         				from being compared to this object.
     */
    public boolean equals(Object obj) {
      return (compareTo((SortContainer) obj) == 0);
    }

    /**
     * Returns a string representation of the sort container.
     *
     * @return		the string representation (value + index)
     */
    public String toString() {
      return "value=" + m_Value + ", index=" + m_Index;
    }
  }
  
  /** the actual table model */
  protected TableModel mModel;

  /** the mapping between displayed and actual index */
  protected int[] mIndices;

  /** the sort column */
  protected int mSortColumn;

  /** whether sorting is ascending or descending */
  protected boolean mAscending;
  
  /**
   * initializes with no model
   */
  public SortedTableModel() {
    this(null);
  }

  /**
   * initializes with the given model
   *
   * @param model       the model to initialize the sorted model with
   */
  public SortedTableModel(TableModel model) {
    setModel(model);
  }

  /**
   * sets the model to use
   *
   * @param value       the model to use
   */
  public void setModel(TableModel value) {
    mModel = value;

    // initialize indices
    if (mModel == null) {
      mIndices = null;
    }
    else {
      initializeIndices();
      mSortColumn = -1;
      mAscending  = true;
      mModel.addTableModelListener(this);
    }
  }

  /**
   * (re-)initializes the indices
   */
  protected void initializeIndices() {
    int       i;

    mIndices = new int[mModel.getRowCount()];
    for (i = 0; i < mIndices.length; i++)
      mIndices[i] = i;
  }

  /**
   * returns the current model, can be null
   *
   * @return            the current model
   */
  public TableModel getModel() {
    return mModel;
  }

  /**
   * returns whether the table was sorted
   *
   * @return        true if the table was sorted
   */
  public boolean isSorted() {
    return (mSortColumn > -1);
  }

  /**
   * whether the model is initialized
   *
   * @return            true if the model is not null and the sort indices
   *                    match the number of rows
   */
  protected boolean isInitialized() {
    return (getModel() != null);
  }

  /**
   * Returns the actual underlying row the given visible one represents. Useful
   * for retrieving "non-visual" data that is also stored in a TableModel.
   * 
   * @param visibleRow	the displayed row to retrieve the original row for
   * @return		the original row
   */
  public int getActualRow(int visibleRow) {
    if (!isInitialized())
      return -1;
    else
      return mIndices[visibleRow];
  }
  
  /**
   * Returns the most specific superclass for all the cell values in the
   * column.
   *
   * @param columnIndex     the index of the column
   * @return                the class of the specified column
   */
  public Class getColumnClass(int columnIndex) {
    if (!isInitialized())
      return null;
    else
      return getModel().getColumnClass(columnIndex);
  }

  /**
   * Returns the number of columns in the model
   *
   * @return          the number of columns in the model
   */
  public int getColumnCount() {
    if (!isInitialized())
      return 0;
    else
      return getModel().getColumnCount();
  }

  /**
   * Returns the name of the column at columnIndex
   *
   * @param columnIndex   the column to retrieve the name for
   * @return              the name of the specified column
   */
  public String getColumnName(int columnIndex) {
    if (!isInitialized())
      return null;
    else
      return getModel().getColumnName(columnIndex);
  }

  /**
   * Returns the number of rows in the model.
   *
   * @return              the number of rows in the model
   */
  public int getRowCount() {
    if (!isInitialized())
      return 0;
    else
      return getModel().getRowCount();
  }

  /**
   * Returns the value for the cell at columnIndex and rowIndex.
   *
   * @param rowIndex      the row
   * @param columnIndex   the column
   * @return              the value of the sepcified cell
   */
  public Object getValueAt(int rowIndex, int columnIndex) {
    if (!isInitialized())
      return null;
    else
      return getModel().getValueAt(mIndices[rowIndex], columnIndex);
  }

  /**
   * Returns true if the cell at rowIndex and columnIndex is editable.
   *
   * @param rowIndex      the row
   * @param columnIndex   the column
   * @return              true if the cell is editable
   */
  public boolean isCellEditable(int rowIndex, int columnIndex) {
    if (!isInitialized())
      return false;
    else
      return getModel().isCellEditable(mIndices[rowIndex], columnIndex);
  }

  /**
   * Sets the value in the cell at columnIndex and rowIndex to aValue.
   *
   * @param aValue        the new value of the cell
   * @param rowIndex      the row
   * @param columnIndex   the column
   */
  public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
    if (isInitialized())
      getModel().setValueAt(aValue, mIndices[rowIndex], columnIndex);
  }

  /**
   * sorts the table over the given column (ascending)
   *
   * @param columnIndex     the column to sort over
   */
  public void sort(int columnIndex) {
    sort(columnIndex, true);
  }

  /**
   * sorts the table over the given column, either ascending or descending
   *
   * @param columnIndex     the column to sort over
   * @param ascending       ascending if true, otherwise descending
   */
  public void sort(int columnIndex, boolean ascending) {
    int       			columnType;
    int       			i;
    ArrayList<SortContainer>	sorted;
    SortContainer		cont;
    Object			value;

    // can we sort?
    if (    (!isInitialized())
         || (getModel().getRowCount() != mIndices.length) ) {

      System.out.println(
          this.getClass().getName() + ": Table model not initialized!");

      return;
    }

    // init
    mSortColumn = columnIndex;
    mAscending  = ascending;
    initializeIndices();
    
    // determine the column type: 0=string/other, 1=comparable
    if (ClassDiscovery.hasInterface(Comparable.class, getColumnClass(mSortColumn)))
      columnType = 1;
    else
      columnType = 0;

    // create list for sorting
    sorted = new ArrayList<SortContainer>();
    for (i = 0; i < getRowCount(); i++) {
      value = mModel.getValueAt(mIndices[i], mSortColumn);
      if (columnType == 0)
	cont = new SortContainer((value == null) ? null : value.toString(), mIndices[i]);
      else
	cont = new SortContainer((Comparable) value, mIndices[i]);
      sorted.add(cont);
    }
    Collections.sort(sorted);

    for (i = 0; i < sorted.size(); i++) {
      if (mAscending)
	mIndices[i] = sorted.get(i).getIndex();
      else
	mIndices[i] = sorted.get(sorted.size() - 1 - i).getIndex();
    }

    sorted.clear();
    sorted = null;
  }

  /**
   * This fine grain notification tells listeners the exact range of cells,
   * rows, or columns that changed.
   *
   * @param e       the event
   */
  public void tableChanged(TableModelEvent e) {
    initializeIndices();
    if (isSorted())
      sort(mSortColumn, mAscending);
    
    fireTableChanged(e);
  }
  
  /**
   * Adds a mouselistener to the header: left-click on the header sorts in
   * ascending manner, using shift-left-click in descending manner.
   *
   * @param table       the table to add the listener to
   */
  public void addMouseListenerToHeader(JTable table) {
    final SortedTableModel modelFinal = this;
    final JTable tableFinal = table;
    tableFinal.setColumnSelectionAllowed(false);
    JTableHeader header = tableFinal.getTableHeader();

    if (header != null) {
      MouseAdapter listMouseListener = new MouseAdapter() {
        public void mouseClicked(MouseEvent e) {
          TableColumnModel columnModel = tableFinal.getColumnModel();
          int viewColumn = columnModel.getColumnIndexAtX(e.getX());
          int column = tableFinal.convertColumnIndexToModel(viewColumn);
          if (    e.getButton() == MouseEvent.BUTTON1 
               && e.getClickCount() == 1 
               && !e.isAltDown()
               && column != -1 ) {
            int shiftPressed = e.getModifiers() & InputEvent.SHIFT_MASK;
            boolean ascending = (shiftPressed == 0);
            modelFinal.sort(column, ascending);
          }
        }
      };
      
      header.addMouseListener(listMouseListener);
    }
  }
}
